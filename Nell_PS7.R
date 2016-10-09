# ROUTE-level data

# lm
#   good Type I error control
#   unbiased
#   this is not linear data
#   poor "interpretability"
# binomial glm
#   better for predictions
#   can give more power, if type I isn't a problem
# quasibinomial glm
#   helpful for overdispersed data
# betabinomial ML with LRT
#   more specifically appropriate for data, given what we know
# betabinomial ML with bootstrapped test
#   as above, plus would give better p-values
# glmm with approximate t-tests (from lmerTest)
#   not sure what the random variable would be
# glmm with bootstrapped test
#   same as above, but bootstrapped test is "gold standard"


# STATION-level data

# lm
#   good type I error control IF stations within routes aren't correlated
# binomial glm
#   better for predictions
# lmm with approximate t-tests (from lmerTest)
#   accounts for autocorrelation within routes
#   good type I error control
# glmm with approximate t-tests (from lmerTest)
#   same as above
#   if type I is okay, can give more power
# Bonus: glmm with bootstrapped test (this could take some computer time)
#   type I should be okay
#   may give more power
#   accounts for within-route autocorrelation






# Install required packages if they're not already installed, then load
for (f in c('magrittr', 'dplyr', 'readr', 'tidyr', 'ggplot2', 'parallel', 'lme4', 
            'car', 'GGally', 'boot')) {
    if (!f %in% rownames(installed.packages())) {
        install.packages(f, dependencies = TRUE)
    }
    library(f, character.only = TRUE)
}; rm(f)

RNGkind("L'Ecuyer-CMRG")


# Inverse logit function
inv_logit <- function(x){
    exp(x)/(1 + exp(x))
}
# This is the logit function
logit <- function(x){
    log(x/(1 - x))
}


# I have 4, but like to leave myself one for other tasks...
ncpus <- 3





# Input to "tidy" data frame with WIND_SPEED^0.5
grouse_df <- read_csv(file = "grouse_data_7Sep16.csv",
                      col_types = paste0(c(rep('?', 12), 'c', rep('?', 17)),
                                         collapse = '')) %>%
    mutate(DATE = as.Date(DATE, format = "%d-%b"),
           WINDSPEEDSQR = sqrt(WIND_SPEED)) %>%
    mutate_each(funs(as.factor), STATION, PERIOD) %>%
    gather(species, detected, RUGR, WITU, STGR) %>%
    mutate(
        species = factor(species, levels = c('STGR', 'RUGR', 'WITU'),
                         labels = c('Sharp-tailed Grouse', 'Ruffed Grouse',
                                    'Wild Turkey'))
    )

rg3 <- grouse_df %>%
    filter(PERIOD == 3, species == 'Ruffed Grouse') %>% 
    select(ROUTE, STATION, detected, WIND_SPEED, WINDSPEEDSQR) %>% 
    rename(RUGR = detected)

rg3_route <- rg3 %>% 
    group_by(ROUTE) %>% 
    summarize(
        RUGR = sum(RUGR),
        STATIONS = n(),
        WIND_SPEED = mean(WIND_SPEED),
        WINDSPEEDSQR = mean(WINDSPEEDSQR))





##################################################################################
# PS7: Maximum Likelihood
##################################################################################


# Questions (which are repeated below)


# 1. Is the ML estimator of b0, b1, or theta biased?

# 2. Explain the pattern of correlation among paramaters. Do these make sense to you?

# 3. Compare the two methods for determining whether b1 differs from zero. Which method 
# is correct?

# 4. Compare the bootstrap tests of H0:b1 = 0 to the LRT. Which do you think gives the 
# best test? What would you want to do to validate the LRT? (Bonus points if you go 
# ahead and do this validation!)

# 5. Perform all of the analyses that we have done so far to test whether WINDSPEEDSQR 
# has an effect on detecting RUGR in PERIOD 3. What are the pros and cons for each 
# analysis? I have given some of the code below that you could use for ML estimation 
# of the beta-binomial and bootstrapping the glmms.



# Description of the data from Michael Hardy

# In 2014, we surveyed 117 roadside routes once during each of three consecutive 18-day
# sampling periods from 31 March to 24 May 2014, the period coinciding with Sharp-tailed 
# Grouse lekking activity. Each route included 8 stations spaced at 0.8 km intervals. At 
# each station, one observer spent 4 min watching and listening for all upland game bird 
# species. Following completion of the survey, the observer recorded environmental 
# covariates that could potentially influence detection probability (wind speed, 
# temperature, humidity, barometric pressure, ambient light, sky code, and various 
# sources of auditory interference). Observers also recorded the presence of raptors 
# and corvid species. Each observer surveyed a set of three sites per morning. 
# Observers were rotated among sites and the order that each site was surveyed (first, 
# second, third) was changed during each of the three sampling periods. Surveys began 
# 45 min before sunrise and were completed by 3 hours after sunrise.




# ========================================================================================
# ========================================================================================

#       ML estimation of the beta-binomial distribution

# ========================================================================================
# ========================================================================================

# --------------------------
# --------------------------
#   Functions
# --------------------------
# --------------------------

# ----------
# Simulation function
# ----------

# Simulate data rg3_route with a beta-binomial at the route level
beta_b_sim <- function(size, X, b0, b1, theta) {
    # Inner function to do a single route's simulation
    route_sim <- function(n, x, b0, b1, theta) {
        p <- inv_logit(b0 + b1 * x)
        prob_route <- rbeta(n = 1, shape1 = (p * theta), shape2 = {(1 - p) * theta})
        rbinom(n = 1, size = n, prob = prob_route)
    }
    # Using `mapply` to iterate through the two input vectors
    mapply(route_sim, size, X, 
           MoreArgs = list(b0 = b0, b1 = b1, theta = theta), 
           USE.NAMES = FALSE, SIMPLIFY = TRUE)
}



# ----------
# Probability distribution and likelihood functions
# ----------

# Probability distribution function for a betabinomial distribution from the 
# library "emdbook"
dbetabinom <- function(y, prob, size, theta, log = FALSE) {
	v <- lfactorial(size) - lfactorial(y) - lfactorial(size - y) - 
	    lbeta(theta * (1 - prob), theta * prob) + 
	    lbeta(size - y + theta * (1 - prob), y + theta * prob)
	if (any(y %% 1 != 0)) {
		warning("non-integer x detected; returning zero probability")
		v[which(y %% 1 != 0)] <- -Inf
	}
	if (log == T) {
		return(v)
	} else {
		return(exp(v))
	}
}


# Log-likelihood function for the betabinomial given data Y (vector of successes), 
# Size (vector of number of trials), and independent variable X (WINDSPEEDSQR) in 
# terms of parameters prob and theta
dbetabinom_LLF <- function(parameters, Y, size, X) {
	theta <- parameters['theta']
	b0 <- parameters['b0']
	b1 <- parameters['b1']

	prob <- inv_logit(b0 + b1 * X)	
	LL <- -sum(dbetabinom(y = Y, size = size, prob = prob, theta = theta, log = TRUE))

	return(LL)
}

# Log-likelihood function for the betabinomial given data Y (vector of successes), 
# Size (vector of number of trials), and independent variable X (WINDSPEEDSQR) with 
# b1 = 0 in terms of parameters prob and theta
dbetabinom_LLF0 <- function(parameters, Y, size, X) {
	theta <- parameters['theta']
	b0 <- parameters['b0']
	
	prob <- inv_logit(b0)
	LL <- -sum(dbetabinom(y = Y, size = size, prob = prob, theta = theta, log = TRUE))

	return(LL)
}


# ----------
# Simulate and fit many times
# ----------

# It is assumed that the optim_fn takes options Y, size, and X, and a vector
# to be optimized over
sim_bb_fit <- function(nsims, optim_fn, optim_par, sim_par, size, X, seed = NULL, 
                       method = "BFGS") {
    sim_b0 <- as.numeric(sim_par['b0'])
    sim_b1 <- as.numeric(sim_par['b1'])
    sim_theta <- as.numeric(sim_par['theta'])
    # Inner function to get parameters b0, b1, and theta once
    # This inherits all arguments from parent environment
    # Argument i is only present so it can be used with lapply
    one_sim <- function(i) {
        sim_y <- beta_b_sim(size, X, sim_b0, sim_b1, sim_theta)
        sim_fit <- suppressWarnings(
            optim(fn = optim_fn, par = optim_par, 
                  Y = sim_y, size = size, X = X, method = method))
        sim_mat <- matrix(sim_fit$par, nrow = 1)
        colnames(sim_mat) <- names(sim_fit$par)
        return(as_data_frame(sim_mat))
    }
    if (!is.null(seed)) set.seed(seed)
    sim_list <- lapply(seq(nsims), one_sim)
    sim_df <- bind_rows(sim_list)
    return(sim_df)
}






# --------------------------
# --------------------------
# Fit to real data
# --------------------------
# --------------------------

rg3_route_fit <- optim(fn = dbetabinom_LLF, par = c(theta = 1, b0 = 0, b1 = 0.1), 
                       Y = rg3_route$RUGR, size = rg3_route$STATIONS, 
                       X = rg3_route$WINDSPEEDSQR, method = "BFGS")

# Distribution of the estimators of b0, b1, and theta

b0 <- rg3_route_fit$par['b0'] %>% as.numeric
b1 <- rg3_route_fit$par['b1'] %>% as.numeric
theta <- rg3_route_fit$par['theta'] %>% as.numeric


sim_bb <- sim_bb_fit(1000, dbetabinom_LLF, c(theta = 1, b0 = 0, b1 = 0.1), 
                     c(theta = theta, b0 = b0, b1 = b1), rg3_route$STATIONS, 
                     rg3_route$WINDSPEEDSQR, seed = 1)


sim_bb %>% 
    gather(parameter, value, everything(), factor_key = TRUE) %>% 
    ggplot(aes(value, fill = parameter)) +
    theme_bw() +
    geom_histogram(bins = 50) +
    geom_vline(
        data = data_frame(x = c(theta, b0, b1), 
                          parameter = factor(c('theta', 'b0', 'b1'))),
        aes(xintercept = x),
        linetype = 3) +
    facet_grid(~ parameter, scales = 'free') +
    theme(legend.position = 'none')

# ----------
# 1. Is the ML estimator of b0, b1, or theta biased?
# ----------

# They do not appear biased from my plots, although theta does not appear to have
# a normal distribution. Presumably, this is because it cannot be < 0.







# ----------
# 2. Explain the pattern of correlation among paramaters. Do these make sense to you?
# ----------

sim_bb %>% 
    ggpairs() +
    theme_bw()


# This makes sense to me that b0 and b1 would be negatively correlated.
# When the intercept is lower, it makes sense that a model with a higher slope would 
# better account for points at higher X and Y values.










# Bootstrap testing H0: b1 = 0 from the distribution of the estimator of b1

ggplot(sim_bb, aes(b1)) +
    geom_histogram(bins = 50, fill = 'dodgerblue') +
    geom_vline(linetype = 3, xintercept = 0) +
    theme_bw() +
    geom_text(data = data_frame(
        x = -0.8, y = 55, 
        label = paste('P-value =', 2 * round(mean(sim_bb$b1 > 0), 3))
    ), 
    aes(x = x, y = y, label = label), hjust = 0, vjust = 1, size = 6) +
    xlab('ML estimates of' ~ beta[1])


# Bootstrap testing H0: b1 = 0 under the null hypothesis

rg3_route_fit0 <- optim(fn = dbetabinom_LLF0, par = c(theta = 1, b0 = 0), 
                       Y = rg3_route$RUGR, size = rg3_route$STATIONS, 
                       X = rg3_route$WINDSPEEDSQR, method = "BFGS")

b0_0 <- rg3_route_fit0$par['b0'] %>% as.numeric
b1_true <- b1
theta_0 <- rg3_route_fit0$par['theta'] %>% as.numeric


sim_bb0 <- sim_bb_fit(1000, dbetabinom_LLF, c(theta = 1, b0 = 0, b1 = 0.1), 
                     c(theta = theta_0, b0 = b0_0, b1 = 0), rg3_route$STATIONS, 
                     rg3_route$WINDSPEEDSQR, seed = 2)




ggplot(sim_bb0, aes(b1)) +
    geom_histogram(bins = 50, fill = 'dodgerblue') +
    geom_vline(linetype = 3, xintercept = b1_true) +
    theme_bw() +
    geom_text(data = data_frame(
        x = 0.5, y = 70, 
        label = paste('P-value =', 2 * round(mean(sim_bb0$b1 < b1_true), 3))
    ), 
    aes(x = x, y = y, label = label), hjust = 1, vjust = 1, size = 6) +
    xlab('ML estimates of' ~ beta[1] * ' under ' * H[0])



# ----------
# 3. Compare the two bootstrap methods for determining whether b1 differs from zero. 
# Which method is correct?
# ----------


# For hypothesis testing, technically the second version is the correct one. The first
# is more appropriate for presenting confidence intervals.









# Likelihood ratio test (LRT)
z <- suppressWarnings(
    optim(fn = dbetabinom_LLF, par = c(theta = 1, b0 = 0, b1 = 0.1), 
          Y = rg3_route$RUGR, size = rg3_route$STATIONS, X = rg3_route$WINDSPEEDSQR, 
          method = "BFGS"))

z0 <- suppressWarnings(
    optim(fn = dbetabinom_LLF0, par = c(theta = 1, b0 = 0), Y = rg3_route$RUGR, 
          size = rg3_route$STATIONS, X = rg3_route$WINDSPEEDSQR, method = "BFGS"))

logLik <- -z$value
logLik0 <- -z0$value

pchisq(2*(logLik - logLik0), df = 1, lower.tail = F)

# ----------
# 4. Compare the bootstrap tests of H0:b1 = 0 to the LRT. Which do you think gives the 
# best test? What would you want to do to validate the LRT? (Bonus points if you go 
# ahead and do this validation!)
# ----------

# The bootstrap test is the "gold standard", but LRT appears to agree quite well.
# To validate, you could run simulations with b1 = 0 and look at type I for LRT.














# ========================================================================================
# ========================================================================================

# Comparison among all methods
# Route level data

# ========================================================================================
# ========================================================================================

# 5. Perform all of the analyses that we have done so far to test whether WINDSPEEDSQR 
# has an effect on detecting RUGR in PERIOD 3. What are the pros and cons for each 
# analysis? I have given some of the code below that you could use for ML estimation of 
# the beta-binomial and bootstrapping the glmms.

# ROUTE-level data

# lm
# binomial glm
# quasibinomial glm
# betabinomial ML with LRT
# betabinomial ML with bootstrapped test
# glmm with approximate t-tests (from lmerTest)
# glmm with bootstrapped test

# STATION-level data

# lm
# binomial glm
# lmm with approximate t-tests (from lmerTest)
# glmm with approximate t-tests (from lmerTest)
# Bonus: glmm with bootstrapped test (this could take some computer time)


##################################################################################
# needed code




# Simulate data w with a logitnormal-binomial at the route level
lnb_sim <- function (size, X, b0, b1, sd) {
    # Inner function to do a single route's simulation
    route_sim <- function(n, x, b0, b1, sd) {
        E <- rnorm(n = 1, mean = 0, sd = sd)
        prob_route <- inv_logit(b0 + b1*x + E)
        rbinom(n = 1, size = n, prob = prob_route)
    }
    # Using `mapply` to iterate through the two input vectors
    mapply(route_sim, size, X, 
           MoreArgs = list(b0 = b0, b1 = b1, sd = sd), 
           USE.NAMES = FALSE, SIMPLIFY = TRUE)
}

# Bootstrap b1 ML estimate
boot_test_b1 <- function(optim_fn, optim_fn_0, Y, size, X, nsims = 1000, seed = NULL,
                         optim_par = c(theta = 1, b0 = 0, b1 = 0.1), 
                         optim_par_0 = c(theta = 1, b0 = 0)) {
    fit_0 <- optim(fn = optim_fn_0, par = optim_par_0, Y = Y, size = size,
                   X = X, method = "BFGS")
    b0_0 <- as.numeric(fit_0$par['b0'])
    theta_0 <- as.numeric(fit_0$par['theta'])
    sim_0 <- sim_bb_fit(nsims, optim_fn, optim_par, 
                        c(theta = theta_0, b0 = b0_0, b1 = 0), size, X, seed = seed)
    return(sim_0$b1)
}




# Extract model info for route-level models
model_ext_route <- function(m, models){
    mmat <- matrix(c(
        m,
        ifelse(m == 'betab', as.numeric(models[[m]][[1]]$par['b1']), 
               models[[m]]$coef[2]),
        ifelse(m == 'betab', 
               pchisq(2 * (as.numeric(-models[[m]][[1]]$value) - 
                               as.numeric(-models[[m]][[2]]$value)), 
                      df = 1, lower.tail = FALSE), 
               summary(models[[m]])$coef[2,4])),
        nrow = 1
    )
    colnames(mmat) <- c('method', 'b1', 'P')
    mdf <- as_data_frame(mmat)
    mdf$b1 <- as.numeric(mdf$b1)
    mdf$P <- as.numeric(mdf$P)
    return(mdf)
}


# *Notes:* It's assumed that the sim_fun takes size_vec, then x_vec, then other args
comp_methods_route <- function(size_vec, x_vec, sim_fun, model_ext_fun, ...) {
    
    sim_y <- sim_fun(size_vec, x_vec, ...)
    
    models <- list(
        lm = lm(sim_y ~ x_vec),
        glm = glm(cbind(sim_y, size_vec - sim_y) ~ x_vec, family = "binomial"),
        qglm = glm(cbind(sim_y, size_vec - sim_y) ~ x_vec, family = "quasibinomial"),
        betab = list(
            suppressWarnings(
                optim(fn = dbetabinom_LLF, par = c(theta = 1, b0 = 0, b1 = 0.1),
                      Y = sim_y, size = size_vec, X = x_vec, method = 'BFGS')),
            suppressWarnings(
                optim(fn = dbetabinom_LLF0, par = c(theta = 1, b0 = 0), Y = sim_y,
                      size = size_vec, X = x_vec, method = "BFGS"))
            )
    )

    lapply(names(models), model_ext_fun, models = models)
}


# Comparing Type I error rates and power (takes ~1.25 min)
b1_range <- seq(0, 0.8, 0.2)
nsims <- 1000

set.seed(3)
rej_route <- mclapply(
    1:nsims,
    function(i) {
        lapply(
            b1_range,
            function(b1_i){
                comp_methods_route(rg3_route$STATIONS, rg3_route$WINDSPEEDSQR, 
                                   beta_b_sim, model_ext_route, b0 = b0_0, b1 = b1_i, 
                                   theta = theta_0)
            }) %>% 
            do.call(what = bind_rows, args = .) %>% 
            mutate(b1_true = rep(b1_range, each = 4))
    }, 
    mc.cores = ncpus) %>%
    bind_rows




rej_route %>% 
    group_by(method, b1_true) %>% 
    summarize(rejected = mean(P < 0.05)) %>% 
    ggplot(aes(b1_true, rejected, color = factor(method))) +
    theme_bw() +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = 3)



# Table of Type I error
rej_route %>% 
    filter(b1_true == 0) %>% 
    group_by(method) %>% 
    summarize(rejected = mean(P < 0.05))

# P value for betabinomial ML bootstrap
set.seed(4)
mean(boot_test_b1(dbetabinom_LLF, dbetabinom_LLF0, rg3_route$RUGR, rg3_route$STATIONS,
             rg3_route$WINDSPEEDSQR, seed = 3, nsims = 1000) < b1_true)


# I don't see any major biases (obv lm isn't interpretable here)
rej_route %>% 
    mutate_each(funs(as.factor), method, b1_true) %>% 
    ggplot(aes(b1)) + 
    geom_histogram(bins = 50, fill = 'dodgerblue') +
    theme_bw() +
    geom_vline(aes(xintercept = as.numeric(paste(b1_true))), linetype = 2) +
    facet_grid(b1_true ~ method, scales = 'free')






# Now checking the same thing if the distribution is from a logitnormal-binomial
# distribution
set.seed(5)
rej_route_ln <- mclapply(
    1:nsims, 
    function(i) {
        lapply(
            b1_range, 
            function(b1_i){
                comp_methods_route(rg3_route$STATIONS, rg3_route$WINDSPEEDSQR, 
                                   lnb_sim, model_ext_route, b0 = b0_0, b1 = b1_i, 
                                   sd = sd(rg3_route$RUGR))
            }) %>% 
            do.call(what = bind_rows, args = .) %>% 
            mutate(b1_true = rep(b1_range, each = 4))
    }, 
    mc.cores = ncpus) %>% 
    bind_rows


rej_route_ln %>% 
    group_by(method, b1_true) %>% 
    summarize(rejected = mean(P < 0.05)) %>% 
    ggplot(aes(b1_true, rejected, color = factor(method))) +
    theme_bw() +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = 3)



# Table of Type I error
rej_route_ln %>% 
    filter(b1_true == 0) %>% 
    group_by(method) %>% 
    summarize(rejected = mean(P < 0.05))



# With higher b1_true, there appears to be more bias
rej_route_ln %>% 
    mutate_each(funs(as.factor), method, b1_true) %>% 
    ggplot(aes(b1)) + 
    geom_histogram(bins = 50, fill = 'dodgerblue') +
    theme_bw() +
    geom_vline(aes(xintercept = as.numeric(paste(b1_true))), linetype = 2) +
    facet_grid(b1_true ~ method, scales = 'free')



# ----------
# Assuming the data follows a betabinomial distribution best, my rankings are as follows:
# lm (most powerful test with good type I error control)
# betabinomial ML with bootstrapped test (if needing predictions, I'd use this one)
# betabinomial ML with LRT
# quasibinomial glm
# binomial glm (would not use at all bc of terrible type I error control)

# Assuming the data follows a logitnormal-binomial distribution best, my rankings are 
# as follows:
# lm (most powerful test with good type I error control and no bias)
# I wouldn't use any of these:
# betabinomial ML with bootstrapped test (biased)
# betabinomial ML with LRT (biased)
# quasibinomial glm (biased)
# binomial glm (terrible type I error control)
# ----------





























# ========================================================================================
# ========================================================================================

# Comparison among all methods
# Station level data

# ========================================================================================
# ========================================================================================


# Simulate data w with a logitnormal-binomial at the station level
lnorm_b_sim <- function(groups, X, b0, b1, sd) {
    # Inner function to do a single group's simulation
    group_sim <- function(X_i) {
        n <- length(X_i)
        E <- rnorm(n = 1, mean = 0, sd = sd)
        prob_obs <- inv_logit(b0 + (b1 * X_i) + E)
        rbinom(n = n, size = 1, prob = prob_obs)
    }

    split_list <- split(X, groups)
    sim_list <- lapply(split_list, group_sim)
    return(as.numeric(c(sim_list, recursive = TRUE)))
}




# Extract model info for route-level models
model_ext_station <- function(m, models){
    if (m == 'lmm') {
        P <- summary(models[[m]])$coef[2,5]
    } else if (m == 'glmm_boot') {
        boot_obj = bootMer(models[[m]], 
                           function(m){summary(m)$coef[2,4]}, 
                           nsim = 100, type = 'parametric')
        P <- 2 * min(mean(boot_obj$t < 0), mean(boot_obj$t > 0))
    } else {
        P <- summary(models[[m]])$coef[2,4]
    }
    b1 <- ifelse(
        grepl('mm', m), 
        summary(models[[m]])$coef[2,1], 
        models[[m]]$coef[2]
        )
    
    mmat <- matrix(c(m, b1, P), nrow = 1)
    colnames(mmat) <- c('method', 'b1', 'P')
    mdf <- as_data_frame(mmat)
    mdf$b1 <- as.numeric(mdf$b1)
    mdf$P <- as.numeric(mdf$P)
    return(mdf)
}





# *Notes:* It's assumed that the sim_fun takes groups_vec, then x_vec, then other args
# It's also assumed that groups_vec is the random effect in mixed models
comp_methods_station <- function(groups_vec, x_vec, sim_fun, model_ext_fun, ...) {
    
    sim_y <- sim_fun(groups_vec, x_vec, ...)
    
    models <- list(
        lm = lm(sim_y ~ x_vec),
        glm = glm(sim_y ~ x_vec, family = "binomial"),
        lmm = lmerTest::lmer(sim_y ~ x_vec + (1 | groups_vec), 
                             control = lmerControl(calc.derivs = FALSE)),
        glmm = glmer(sim_y ~ x_vec + (1 | groups_vec), family = "binomial", nAGQ = 0,
                     control = glmerControl(calc.derivs = FALSE))
        )
    models[['glmm_boot']] <- models[['glmm']]
    
    lapply(names(models), model_ext_fun, models = models)
}



# system.time(
# comp_methods_station(rg3$ROUTE, rg3$WINDSPEEDSQR, lnorm_b_sim,
#                       model_ext_station, b0 = mean(rg3_route$RUGR/rg3_route$STATIONS),
#                       b1 = 0,
#                       sd = sd(rg3_route$RUGR))
# )
# # Time taken per iteration
# #    user  system elapsed 
# #   4.875   0.108   5.030

# Comparing Type I error rates and power (takes 3.715058 hours)
b1_range <- seq(0, 0.5, 0.1)
nsims <- 1000


# set.seed(11)
# rej_station <- mclapply(
#     1:nsims,
#     function(i) {
#         lapply(
#             b1_range,
#             function(b1_i){
#                 comp_methods_station(rg3$ROUTE, rg3$WINDSPEEDSQR, lnorm_b_sim,
#                                      model_ext_station, 
#                                      b0 = mean(rg3_route$RUGR/rg3_route$STATIONS),
#                                      b1 = b1_i,
#                                      sd = sd(rg3_route$RUGR))
#             }) %>% 
#             do.call(what = bind_rows, args = .) %>% 
#             mutate(b1_true = rep(b1_range, each = 5))
#     }, 
#     mc.cores = ncpus) %>%
#     bind_rows
# 
# save(rej_station, file = 'station_sims.RData', compression_level = 9)
load('station_sims.RData')


rej_station %>% 
    group_by(method, b1_true) %>% 
    summarize(rejected = mean(P < 0.05)) %>% 
    ggplot(aes(b1_true, rejected, color = factor(method), linetype = factor(method))) +
    theme_bw() +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = 3)



# Table of Type I error
rej_station %>% 
    filter(b1_true == 0) %>% 
    group_by(method) %>% 
    summarize(rejected = mean(P < 0.05))



# 
# # Moving this file to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp Nell_PS7.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/'")
# )
