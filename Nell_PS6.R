
# Install required packages if they're not already installed, then load
for (f in c('magrittr', 'dplyr', 'readr', 'tidyr', 'ggplot2', 'parallel', 'lme4', 
            'car')) {
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


# I have 4, but like to leave myself one for other tasks...
ncpus <- 3




# 
# # Input to "tidy" data frame with WIND_SPEED^0.5
# grouse_df <- read_csv(file = "grouse_data_7Sep16.csv", 
#                       col_types = paste0(c(rep('?', 12), 'c', rep('?', 17)), 
#                                          collapse = '')) %>%
#     mutate(DATE = as.Date(DATE, format = "%d-%b"),
#            WINDSPEEDSQR = sqrt(WIND_SPEED)) %>%
#     mutate_each(funs(as.factor), STATION, PERIOD) %>%
#     gather(species, detected, RUGR, WITU, STGR) %>%
#     mutate(
#         species = factor(species, levels = c('STGR', 'RUGR', 'WITU'), 
#                          labels = c('Sharp-tailed Grouse', 'Ruffed Grouse', 
#                                     'Wild Turkey'))
#     )
# 
# 
# g3 <- grouse_df %>% 
#     filter(PERIOD == 3)






##################################################################################
# PS6: LMMs and GLMMs
##################################################################################
# Questions (which are repeated below):

# 1. Compare three methods for analyzing a linear model with 10 routes and 20 stations 
# per route, with routes having the same intercept but different slopes: the methods 
# are (i) 10 separate linear models, (ii) a single linear model with categorical terms 
# for route and route:slope interactions, and (iii) a linear mixed model with route as 
# a random effect for both slope and intercept. For the simulations, set b0.mean <- 1, 
# b0.sd <- 0, b1.mean <- .2, b1.sd <- .05, and e.sd <- 0.1. What are the mean and 
# variance in slopes calculated with each of the three methods?

# 2. Compare two methods for analyzing a linear model with 10 routes and 20 stations 
# per route, with routes having the same slope but different intercepts: the methods 
# are (i) a single linear model with categorical terms for route, and (ii) a linear 
# mixed model with route as a random effect. For the simulations, set b0.mean <- 1, 
# b0.sd <- 1, b1.mean <- .2, b1.sd <- 0, and e.sd <- 1. Note that you will also have 
# to change the models for fitting the data. Which method gives on average the lowest 
# P-values for the significance of the slope? You will have to run the simulations 
# several times for this.

# 3. Compare three methods for analyzing a binary model with 10 routes and 20 stations 
# per route, with routes having different slopes: the methods are (i) 10 separate 
# linear models, (ii) a single linear model with categorical terms for route and 
# route:slope interactions, and (iii) a linear mixed model with route as a random 
# effect for both slope and intercept. For the simulations, set b0.mean <- 0, 
# b0.sd <- 0, b1.mean <- 1, and b1.sd <- .5. What are the mean and variance in slopes 
# calculated with each of the three methods? Are the patterns you get different from 
# those in questions #1 for linear models?

# 4. Compare two methods for analyzing a binary model with 10 routes and 20 stations per 
# route, with routes having the same slope but different intercepts: the methods are (i) 
# a single linear model with categorical terms for route, and (ii) a linear mixed model 
# with route as a random effect. For the simulations, set b0.mean <- 0, b0.sd <- .5, 
# b1.mean <- .2, and b1.sd <- 0. Note that you will also have to change the models for 
# fitting the data. Which method gives on average the lowest P-values for the 
# significance of the slope? You will have to run the simulations several times for this.

# 5. Bonus (optional) question: Write code and check the type I error rate for the glmm 
# testing the hypothesis that there is no variance in the slopes among routes when 
# there is variation in the intercept.


####################################################################
# NOTE: lmer() likes it better with nstations=20, so I changed it to this. But there are
# still sometimes error messages. Just ignore them.

# LMM with variation in the intercept and slope


# =====================================
# =====================================
# 
# 1. Compare three methods for analyzing a linear model with 10 routes and 20 stations 
# per route, with routes having the same intercept but different slopes: the methods 
# are (i) 10 separate linear models, (ii) a single linear model with categorical terms 
# for route and route:slope interactions, and (iii) a linear mixed model with route as 
# a random effect for both slope and intercept. For the simulations, set b0.mean <- 1, 
# b0.sd <- 0, b1.mean <- .2, b1.sd <- .05, and e.sd <- 0.1. What are the mean and 
# variance in slopes calculated with each of the three methods?
# 
# =====================================
# =====================================

norm_sim <- function(i, nstations = 20, 
                    b0.mean = 1, b0.sd = 0, b1.mean = 0.2, b1.sd = 0.05, 
                    e.sd = 0.1, x.mean = 0, x.sd = 1) {
    b0 <- b0.mean + rnorm(n=1, mean=0, sd=b0.sd)
    b1 <- b1.mean + rnorm(n=1, mean=0, sd=b1.sd)
    X <- rnorm(n=nstations, mean=x.mean, sd=x.sd)
    Y <- b0 + b1 * X + rnorm(n=nstations, mean=0, sd=e.sd)
    return(data_frame(X, Y, route = i))
}
mult_norm_sims <- function(nsims, seed = NULL, nroutes = 10, nstations = 20, ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    sim_list <- replicate(nsims, lapply(seq(nroutes), norm_sim, nstations, ...), 
                          simplify = FALSE)
    sim_df <- do.call(bind_rows, sim_list)
    sim_df$sim <- rep(seq(nsims), each = {nstations * nroutes})
    return(sim_df)
}


nroutes <- 10
d <- mult_norm_sims(1, 999)

# -------------
# (i) 10 separate linear models
# -------------
b0.est.list <- 0
b1.est.list <- 0
for(i in 1:nroutes) {
	z <- lm(Y ~ X, data = d[d$route == i,])
	b0.est.list[i] <- z$coef[1]
	b1.est.list[i] <- z$coef[2]
}

# Mean and sd of coefficients
c(mean(b0.est.list), sd(b0.est.list))
c(mean(b1.est.list), sd(b1.est.list))

# -------------
# (ii) a single linear model with categorical terms for route and route:slope interactions
# -------------
z.lm <- lm(Y ~ 0 + X + factor(route) + X:factor(route), data = d)
# summary(z.lm)
# Anova(z.lm)

# Mean and sd of b1
c(mean(z.lm$coef[1] + c(0, z.lm$coef[12:20])), sd(c(0, z.lm$coef[12:20])))


# -------------
# (iii) a linear mixed model with route as a random effect for both slope and intercept
# -------------
z.lmm <- lmer(Y ~ X + (1 + X | route), data = d)
summary(z.lmm)
summary(z.lmm)$coefficients

coef(z.lmm)$route$X %T>% {print(mean(.)); .} %>% sd


# Test the significance of the slope random effect
z0.lmm <- lmer(Y ~ X + (1 | route), data = d)
P.value <- pchisq(2*(logLik(z.lmm) - logLik(z0.lmm)), df=2, lower.tail = F)
P.value


# -------------
# 1. Compare three methods for analyzing a linear model ...
# with routes having the same intercept but different slopes... 
# What are the mean and variance in slopes calculated with each of the three methods?
# -------------

# (i) 10 separate linear models
# > b1.est.list %T>% {print(mean(.)); .} %>% sd
# [1] 0.1748774
# [1] 0.07008706


# (ii) a single linear model with categorical terms for route and route:slope interactions
# > c(mean(z.lm$coef[1] + c(0, z.lm$coef[12:20])), sd(c(0, z.lm$coef[12:20])))
# [1] 0.17487744 0.07008706


# (iii) a linear mixed model with route as a random effect for both slope and intercept
# > coef(z.lmm)$route$X %T>% {print(mean(.)); .} %>% sd
# [1] 0.1760874
# [1] 0.05883794



# =====================================
# =====================================
# 
# 2. Compare two methods for analyzing a linear model with 10 routes and 20 stations 
# per route, with routes having the same slope but different intercepts: the methods 
# are (i) a single linear model with categorical terms for route, and (ii) a linear 
# mixed model with route as a random effect. For the simulations, set b0.mean <- 1, 
# b0.sd <- 1, b1.mean <- .2, b1.sd <- 0, and e.sd <- 1. Note that you will also have 
# to change the models for fitting the data. Which method gives on average the lowest 
# P-values for the significance of the slope? You will have to run the simulations 
# several times for this.
# 
# =====================================
# =====================================



dd <- mult_norm_sims(10, 888, b0.mean = 1, b0.sd = 1, b1.mean = 0.2, b1.sd = 0, e.sd = 1)

# -------------
# (i) a single linear model with categorical terms for route
# -------------
lm_get_P <- function(dat) {
    lm_mod <- lm(Y ~ 0 + X + factor(route), data = dat)
    P <- Anova(lm_mod)$`Pr(>F)`[1]
    return(data_frame(P))
}

lm_P <- dd %>% 
    group_by(sim) %>% 
    do(lm_get_P(.))


# -------------
# (ii) a linear mixed model with route as a random effect.
# -------------

lmm_get_P <- function(dat) {
    lmm_mod <- lmer(Y ~ X + (1 | route), data = dat)
    P <- Anova(lmm_mod)$`Pr(>Chisq)`
    return(data_frame(P))
}


lmm_P <- dd %>% 
    group_by(sim) %>% 
    do(lmm_get_P(.))


# -------------
# 2. Compare two methods for analyzing a linear model ...
# with routes having the same slope but different intercepts
# Which method gives on average the lowest P-values for the significance of the slope?
# -------------


# mean(lm_P$P); mean(lmm_P$P)
# [1] 3.668048e-56
# [1] 4.800958e-116

# The LMM has lower P values














####################################################################
# Binary GLMM with variation in the intercept and slope


# =====================================
# =====================================
# 
# 3. Compare three methods for analyzing a binary model with 10 routes and 20 stations 
# per route, with routes having different slopes: the methods are (i) 10 separate 
# linear models, (ii) a single linear model with categorical terms for route and 
# route:slope interactions, and (iii) a linear mixed model with route as a random 
# effect for both slope and intercept. For the simulations, set b0.mean <- 0, 
# b0.sd <- 0, b1.mean <- 1, and b1.sd <- .5. What are the mean and variance in slopes 
# calculated with each of the three methods? Are the patterns you get different from 
# those in questions #1 for linear models?
# 
# =====================================
# =====================================

nroutes <- 10

# Simulation of binary data
bin_sim <- function(i, nstations = 20, 
                    b0.mean = 0, b0.sd = 0, b1.mean = 1, b1.sd = 0.5, 
                    e.sd = 0.1, x.mean = 0, x.sd = 1) {
    b0 <- b0.mean + rnorm(n=1, mean=0, sd=b0.sd)
    b1 <- b1.mean + rnorm(n=1, mean=0, sd=b1.sd)
    X <- rnorm(n=nstations, mean=x.mean, sd=x.sd)
    p <- inv_logit(b0 + b1 * X)
    Y <- rbinom(n = nstations, size = 1, prob = p)
    return(data_frame(X, Y, route = i))
}

mult_bins_sims <- function(nsims, seed = NULL, nroutes = 10, nstations = 20, ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    sim_list <- replicate(nsims, lapply(seq(nroutes), bin_sim, nstations, ...), 
                          simplify = FALSE)
    sim_df <- do.call(bind_rows, sim_list)
    sim_df$sim <- rep(seq(nsims), each = {nstations * nroutes})
    return(sim_df)
}


d <- mult_bins_sims(1, 999)



# -------------
# (i) 10 separate linear models
# -------------
b0.est.list <- 0
b1.est.list <- 0
for(i in 1:nroutes) {
	z <- glm(Y ~ X, family = "binomial", data = d[d$route == i,])
	b0.est.list[i] <- z$coef[1]
	b1.est.list[i] <- z$coef[2]
}

# Calculate the mean and sd of b1 from the separate GLMs
b1.mean.sd <- c(mean(b1.est.list), sd(b1.est.list))
b1.mean.sd

# -------------
# (ii) a single linear model with categorical terms for route and route:slope interactions
# -------------
# Estimate GLM with a categorical variable for route affecting the slope and the intercept
z.glm <- glm(Y ~ 0 + X + factor(route) + X:factor(route), family = "binomial", data = d)

# b0.mean.sd <- c(mean(z.glm$coef[2:11]), sd(z.glm$coef[2:11]))
# b0.mean.sd
b1.mean.sd <- c(mean(z.glm$coef[1] + c(0, z.glm$coef[12:20])), 
                sd(c(0, z.glm$coef[12:20])))
b1.mean.sd

# -------------
# (iii) a linear mixed model with route as a random effect for both slope and intercept
# -------------
# Estimate GLMM with route as a random effect for the slope and the intercept
z.glmm <- glmer(Y ~ X + (1 + X | route), family = "binomial", data = d, 
                glmerControl(calc.derivs = F))
# summary(z.glmm)

c(mean(coef(z.glmm)$route[,2]), sd(coef(z.glmm)$route[,2]))


# Test the significance of the slope random effect
z0.glmm <- glmer(Y ~ X + (1 | route), family = "binomial", data = d, 
                 glmerControl(calc.derivs=F))
P.value <- pchisq(2*(as.numeric(logLik(z.glmm)) - as.numeric(logLik(z0.glmm))), 
                  df=2, lower.tail = F)
P.value



# -------------
# 3. Compare three methods for analyzing a binary model ...
# with routes having different slopes... 
# What are the mean and variance in slopes calculated with each of the three methods? 
# Are the patterns you get different from those in questions #1 for linear models?
# -------------

# For (i) 10 separate linear models:
# > b1.mean.sd
# [1] 0.8878563 0.5493154


# For (ii), a glm with categorical variable for route and route:slope interactions:
# > b1.mean.sd
# [1] 0.8878563 0.5493154


# For (iii) a glmm with route as a random effect for both slope and intercept:
# c(mean(coef(z.glmm)$route[,2]), sd(coef(z.glmm)$route[,2]))
# [1] 0.7533037 0.0000000


# As in the linear model, the variance for the mixed model version is lower than either
# non-mixed model. The actual mean and sd are actually more incorrect for the mixed model,
# which is perplexing...





# =====================================
# =====================================
# 
# 4. Compare two methods for analyzing a binary model with 10 routes and 20 stations per 
# route, with routes having the same slope but different intercepts: the methods are (i) 
# a single linear model with categorical terms for route, and (ii) a linear mixed model 
# with route as a random effect. For the simulations, set b0.mean <- 0, b0.sd <- .5, 
# b1.mean <- .2, and b1.sd <- 0. Note that you will also have to change the models for 
# fitting the data. Which method gives on average the lowest P-values for the 
# significance of the slope? You will have to run the simulations several times for this.
# 
# =====================================
# =====================================


dd <- mult_bins_sims(10, 888, b0.mean = 0, b0.sd = 0.5, b1.mean = 0.2, b1.sd = 0)

# -------------
# (i) a single linear model with categorical terms for route
# -------------
glm_get_P <- function(dat) {
    glm_mod <- glm(Y ~ 0 + X + factor(route), family = "binomial", data = dat)
    P <- Anova(glm_mod, test.statistic = 'F')$`Pr(>F)`[1]
    return(data_frame(P))
}

glm_P <- dd %>% 
    group_by(sim) %>% 
    do(glm_get_P(.))




# -------------
# (ii) a linear mixed model with route as a random effect.
# -------------

glmm_get_P <- function(dat) {
    glm_mod <- glmer(Y ~ 0 + X + (1 | route), family = "binomial", data = dat, 
                     glmerControl(calc.derivs = FALSE))
    P <- Anova(glm_mod)$`Pr(>Chisq)`
    return(data_frame(P))
}

glmm_P <- dd %>% 
    group_by(sim) %>% 
    do(glmm_get_P(.))


# -------------
# 4. Compare two methods for analyzing a linear model ...
# with routes having the same slope but different intercepts
# Which method gives on average the lowest P-values for the significance of the slope?
# -------------


# mean(glm_P$P); mean(glmm_P$P)
# [1] 0.0001272579
# [1] 0.0001814222

# The GLMM has, on average, slightly lower P values








# =====================================
# =====================================
# 
# 5. Bonus (optional) question: Write code and check the type I error rate for the glmm 
# testing the hypothesis that there is no variance in the slopes among routes, when 
# there is variation in the intercept.
# 
# =====================================
# =====================================


dd <- mult_bins_sims(1000, 7, b0.mean = 0, b0.sd = 1, b1.mean = 0, b1.sd = 0)

glmm_get_rP <- function(dat) {
    z.sim <- glmer(Y ~ X + (1 + X | route), family="binomial", 
                   data = dat, control = glmerControl(calc.derivs = FALSE))
    z1.sim <- glmer(Y ~ X + (1 | route) + (0 + X | route), family = "binomial", 
                    data = dat, control = glmerControl(calc.derivs = FALSE))
    z0.sim <- glmer(Y ~ X + (1 | route), family = "binomial", 
                    data = dat, control = glmerControl(calc.derivs = FALSE))
    P.b1.sd <- 0.5 * pchisq(2*(logLik(z.sim) - logLik(z0.sim)), df=1, lower.tail = F) + 
        0.5 * pchisq(2*(logLik(z.sim) - logLik(z0.sim)), df=2, lower.tail = F)
    P.b1.sd1 <- .5*pchisq(2*(logLik(z1.sim) - logLik(z0.sim)), df=1, lower.tail = F)
    return(data_frame(P.b1.sd, P.b1.sd1))
}

glmm_rP <- dd %>% 
    group_by(sim) %>% 
    do(glmm_get_rP(.))

mean(glmm_rP$P.b1.sd < 0.05)
# [1] 0.179
mean(glmm_rP$P.b1.sd1 < 0.05)
# [1] 0.264



# 
# # Moving this file to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp Nell_PS6.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/'")
# )
