
# Install required packages if they're not already installed, then load
for (f in c('magrittr', 'dplyr', 'readr', 'tidyr', 'ggplot2', 'parallel')) {
    if (!f %in% rownames(installed.packages())) {
        install.packages(f, dependencies = TRUE)
    }
    library(f, character.only = TRUE)
}; rm(f)


# Inverse logit function
inv_logit <- function(x){
    exp(x)/(1 + exp(x))
}
# This is the logit function
logit <- function(x){
    log(x/(1 - x))
}




##################################################################################
# PS4: Simulation
##################################################################################

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

#           Read and visualize data

# ========================================================================================
# ========================================================================================

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


g3 <- grouse_df %>% 
    filter(PERIOD == 3)


ggplot(g3, aes(WINDSPEEDSQR)) + geom_histogram(bins = 20)

g3 %>% 
    filter(species == 'Ruffed Grouse') %T>% 
    {plot(ggplot(data = ., aes(WINDSPEEDSQR, detected)) + geom_point()); .} %>% 
    ggplot(aes(factor(detected), WINDSPEEDSQR)) + geom_boxplot() + coord_flip()



g3_by_route <- g3 %>% 
    filter(species == 'Ruffed Grouse') %>% 
    select(detected, ROUTE, STATION, WIND_SPEED, WINDSPEEDSQR) %>% 
    group_by(ROUTE) %>% 
    summarize(
        RUGR = sum(detected),
        STATIONS = n(),
        WIND_SPEED = mean(WIND_SPEED),
        WINDSPEEDSQR = mean(WINDSPEEDSQR))

g3_by_route %>% 
    ggplot(aes(WINDSPEEDSQR, RUGR/STATIONS)) + 
    geom_point() + 
    ggtitle('Data aggregated to ROUTES')




# ========================================================================================
# ========================================================================================

#           Beta binomial

# ========================================================================================
# ========================================================================================


# Simulate data w with a beta-binomial at the route level
beta_b_sim <- function(ns, xs, b0, b1, theta) {
    # Inner function to do a single route's simulation
    route_sim <- function(n, x, b0, b1, theta) {
        p <- inv_logit(b0 + b1 * x)
        prob_route <- rbeta(n = 1, shape1 = (p * theta), shape2 = {(1 - p) * theta})
        rbinom(n = 1, size = n, prob = prob_route)
    }
    # Using `mapply` to iterate through the two input vectors
    mapply(route_sim, ns, xs, 
           MoreArgs = list(b0 = b0, b1 = b1, theta = theta), 
           USE.NAMES = FALSE, SIMPLIFY = TRUE)
}



# -----------------------
# 1. What values of p and theta of the beta distribution give simulations most similar
# to the data without including wind speed (i.e., b1 = 0)?
# -----------------------



p_thetas <- expand.grid(theta = c(10, 25, 50, 100), p = seq(0.05, 0.25, 0.05))

set.seed(111)
pt_sims <- mapply(function(p, t) {
        beta_b_sim(g3_by_route$STATIONS, g3_by_route$WINDSPEEDSQR, logit(p), 0, t)
    }, p_thetas$p, p_thetas$theta) %>% 
    as_data_frame() %>% 
    mutate(STATIONS = g3_by_route$STATIONS) %>% 
    gather(SCENARIO, RUGR, starts_with('V')) %>%
    mutate(SCENARIO = sapply(SCENARIO, 
                             function(x) {
                                 i <- as.integer(gsub('V', '', x))
                                 paste('p =', p_thetas$p[i], 
                                       ' theta =', p_thetas$theta[i])
                             }) %>% 
               factor(., levels = paste('p =', p_thetas$p, ' theta =', p_thetas$theta)))




pt_sims %>% 
    ggplot(aes(RUGR/STATIONS)) +
    geom_density(data = g3_by_route, color= 'red') +
    geom_density() +
    theme_bw() +
    facet_wrap(~SCENARIO, ncol = 4, scales = 'free_y')


# From these, I'd say p = 0.15 and theta = 10, 
# although p = 0.2 and theta = 10 looks pretty good, too







# ========================================================================================
# ========================================================================================

#           Logitnormal-Binomial

# ========================================================================================
# ========================================================================================


# Simulate data w with a logitnormal-binomial at the station level
lnorm_b_sim <- function(groups, xs, b0, b1, sd) {
    # Inner function to do a single group's simulation
    group_sim <- function(group_i) {
        xs_i <- xs[groups == group_i]
        n <- length(xs_i)
        prob_route <- rnorm(n = 1, mean = 0, sd = sd)
        prob_obs <- inv_logit(b0 + (b1 * xs_i) + prob_route)
        rbinom(n = n, size = 1, prob = prob_obs)
    }
    
    if (class(groups) == 'factor') {
        unq_groups <- levels(groups)
    } else {
        unq_groups <- unique(groups)
    }
    
    c(sapply(unq_groups, group_sim, USE.NAMES = FALSE), recursive = TRUE)
}


# -----------------------
# 2. What values of b0 and sd.RUGR of the logitnormal distribution give simulations most
# similar to the data without including wind speed (i.e., b1 = 0)?
# -----------------------

b0_sds <- expand.grid(b0 = seq(-2, -10, -2), sd = seq(5))


set.seed(222)
bs_sims <- mapply(function(b, s) { lnorm_b_sim(g3$ROUTE, g3$WINDSPEEDSQR, b, 0, s) },
                  b0_sds$b0, b0_sds$sd) %>% 
    as_data_frame() %>% 
    mutate(ROUTE = g3$ROUTE) %>% 
    gather(SCENARIO, RUGR, starts_with('V')) %>%
    mutate(SCENARIO = sapply(SCENARIO, 
                             function(x) {
                                 i <- as.integer(gsub('V', '', x))
                                 paste('b0 =', b0_sds$b0[i], 
                                       ' sd =', b0_sds$sd[i])
                             }) %>% 
               factor(., levels = paste('b0 =', b0_sds$b0, ' sd =', b0_sds$sd)))


bs_sims %>% 
    group_by(SCENARIO, ROUTE) %>% 
    summarize(
        RUGR = sum(RUGR),
        STATIONS = n()) %>% 
    ungroup %>% 
    ggplot(aes(RUGR/STATIONS)) +
    geom_density(data = g3_by_route, color= 'red') +
    geom_density() +
    theme_bw() +
    facet_wrap(~SCENARIO, ncol = 5, scales = 'free_y')



# b0 = -4, sd = 4 seems to be the best fit
















# ========================================================================================
# ========================================================================================

# Estimating models fit to the simulation data sets

# ========================================================================================
# ========================================================================================


# Function to compare GLM, LM, and quasiGLM for a simulation

comp_methods <- function(x_vec, n_vec, sim_fun, ...) {
    
    sim_y <- sim_fun(...)
    
    models <- list(
        lm = lm(sim_y ~ x_vec),
        glm = glm(cbind(sim_y, n_vec - sim_y) ~ x_vec, family = "binomial"),
        glm.quasi = glm(cbind(sim_y, n_vec - sim_y) ~ x_vec, family = "quasibinomial")
    )
    
    lapply(names(models), 
           function(m){
               data_frame(
                   method = m,
                   b0.est = models[[m]]$coef[1],
                   b1.est = models[[m]]$coef[2],
                   P = summary(models[[m]])$coef[2,4]
               )
           }
    ) %>% bind_rows
}



# ========================================================================================
# ========================================================================================

#    Beta binomial

# ========================================================================================
# ========================================================================================


# -----------------------
# 3. Investigate the distribution of the LM, binomial GLM, and quasibinomial GLM 
# estimators of b1. Do they estimate the same thing? Are the GLM estimators biased? 
# What happens tot he distribution of the estimators as you increase the 
# route-to-route variation (decrease theta)? Note: for b0 use the value that best 
# fits the data (Q1)
# -----------------------




p <- 0.15
theta <- 10
b1 <- 0


nsims <- 1000

set.seed(333)
comp_beta_b <- replicate(nsims, 
                         comp_methods(g3_by_route$WINDSPEEDSQR, g3_by_route$STATIONS,
                                      beta_b_sim, 
                                      xs = g3_by_route$WINDSPEEDSQR, 
                                      ns = g3_by_route$STATIONS, 
                                      b0 = logit(p), b1 = b1, theta = theta), 
                         simplify = FALSE) %>% 
    bind_rows %>% 
    mutate(method = as.factor(method))

comp_beta_b %>% 
    ggplot(aes(b1.est, color = method, linetype = method)) +
    geom_freqpoly(bins = 50, size = 0.75) +
    geom_vline(xintercept = 0, linetype = 3) +
    theme_bw() +
    scale_linetype_manual(values = c(1, 3, 2))



# They appear to estimate the same thing, and I don't detect any major biases


# Changing theta
set.seed(444)
theta_sims <- lapply(
    seq(0, 10, 5), 
    function(t){
        replicate(nsims, 
                  comp_methods(g3_by_route$WINDSPEEDSQR, g3_by_route$STATIONS,
                               beta_b_sim, 
                               xs = g3_by_route$WINDSPEEDSQR, 
                               ns = g3_by_route$STATIONS, 
                               b0 = logit(p), b1 = b1, theta = t), 
                  simplify = FALSE) %>% 
            bind_rows %>% 
            mutate(theta = t)
    }
    ) %>% 
    bind_rows %>% 
    mutate_each(funs(as.factor), method, theta)



theta_sims %>% 
    ggplot(aes(b1.est, color = method, linetype = method)) +
    geom_freqpoly(bins = 50, size = 0.75) +
    geom_vline(xintercept = 0, linetype = 3) +
    theme_bw() +
    scale_linetype_manual(values = c(1, 3, 2)) +
    facet_wrap( ~ theta, ncol = 1, scales = 'free_y')

# At theta = 0, both glm models don't have as much variance in b1.est as the linear model




# -----------------------
# 4. Check the type I errors. Does type I error control for LM, binomial GLM, and 
# quasibinomial GLM estimators depend on theta? How would this affect your 
# interpretation of p-values in your analyses?
# -----------------------


theta_sims %>% 
    group_by(theta, method) %>% 
    summarize(typeI = mean(P < 0.05)) %>% 
    ungroup %>% 
    mutate(theta = as.integer(paste(theta))) %T>% 
    print(.) %>% 
    ggplot(aes(theta, typeI, color = method)) + 
    geom_line(size = 1) +
    theme_bw() +
    geom_hline(yintercept = 0.05, linetype = 3)


# GLM has way higher type I error than the other two, and that increases with lower theta

# You'd have to lower your cutoff if using a GLM:

theta_sims %>% 
    filter(method == 'glm') %>% 
    group_by(theta) %>% 
    summarize(alpha = quantile(P, probs = 0.05))






# -----------------------
# 5. Compare the power of GLMs and LMs with the same data sets (i.e., in the loop over
# simulated data sets, test the hypothesis H0: b1 = 0). Which estimator is more powerful? 
# -----------------------

p <- 0.15
theta <- 10


# Changing b1
set.seed(555)
b1_sims <- lapply(
    seq(0, 0.5, 0.1), 
    function(b){
        replicate(nsims, 
                  comp_methods(g3_by_route$WINDSPEEDSQR, g3_by_route$STATIONS,
                               beta_b_sim, 
                               xs = g3_by_route$WINDSPEEDSQR, 
                               ns = g3_by_route$STATIONS, 
                               b0 = logit(p), b1 = b, theta = theta), 
                  simplify = FALSE) %>% 
            bind_rows %>% 
            mutate(b1 = b)
    }
    ) %>% 
    bind_rows %>% 
    mutate_each(funs(as.factor), method, b1)


b1_sims %>% 
    group_by(b1, method) %>% 
    summarize(rejected = mean(P < 0.05)) %>% 
    ungroup %>% 
    mutate(b1 = as.numeric(paste(b1))) %T>% 
    print(.) %>% 
    ggplot(aes(b1, rejected, color = method)) + 
    geom_line(size = 1) +
    theme_bw()


# glm is technically more powerful from this output, but if you include the type I error
# problems the glm had, then the lm is most powerful









# ========================================================================================
# ========================================================================================

# Logitnormal-binomial

# ========================================================================================
# ========================================================================================



# -----------------------
# 6. Investigate the distribution of the LM, binomial GLM, and quasibinomial GLM 
# estimators of b1. Do they estimate the same thing? Are the GLM estimators biased? 
# What happens as you increase the route-to-route variation (increase sd.RUGR)? 
# Note: for b0 use the value that best fits the data (Q2)
# -----------------------

b0_true <- -4
sd_true <- 4

nsims <- 1000



set.seed(666)
sd_sims <- lapply(
    seq(sd_true, sd_true + (5*4), 5), 
    function(sd_i){
        replicate(nsims, 
                  comp_methods(g3$WINDSPEEDSQR, rep(1, nrow(g3)),
                               lnorm_b_sim, 
                               groups = g3$ROUTE, 
                               xs = g3$WINDSPEEDSQR, 
                               b0 = b0_true, b1 = 0, sd = sd_i), 
                  simplify = FALSE) %>% 
            bind_rows %>% 
            mutate(sd.true = sd_i)
    }
    ) %>% 
    bind_rows %>% 
    mutate_each(funs(as.factor), method, sd.true)


sd_sims %>% 
    ggplot(aes(b1.est, color = method, linetype = method)) +
    geom_freqpoly(bins = 50, size = 0.75) +
    geom_vline(xintercept = 0, linetype = 3) +
    theme_bw() +
    scale_linetype_manual(values = c(1, 3, 2)) +
    facet_wrap( ~ sd.true, ncol = 1, scales = 'free_y')



# For all values of sd, the variance in b1 estimates is way higher for glms than
# the lm





# -----------------------
# 7. Check the type I errors. Does type I error control for LM and GLM estimators 
# depend on sd.RUGR? How would this affect your interpretation of p-values in your 
# analyses?
# -----------------------



sd_sims %>% 
    group_by(sd.true, method) %>% 
    summarize(typeI = mean(P < 0.05)) %>% 
    ungroup %>% 
    mutate(sd.true = as.numeric(paste(sd.true))) %T>% 
    print(.) %>% 
    ggplot(aes(sd.true, typeI, color = method)) + 
    geom_line(size = 1) +
    theme_bw() +
    geom_hline(yintercept = 0.05, linetype = 3)


# All models have *way* higher type I error than they should, and that increases with 
# higher sd

# You'd have to lower your cutoff if using this distribution:

sd_sims %>% 
    filter(method == 'glm') %>% 
    group_by(sd.true) %>% 
    summarize(alpha = quantile(P, probs = 0.05))



# This is so bad that it makes me think I coded this wrong...







# # Moving this file to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp Nell_PS4.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/'")
# )

