
# Install required packages if they're not already installed, then load
for (f in c('magrittr', 'dplyr', 'readr', 'tidyr', 'ggplot2', 'parallel', 'lme4', 
            'lmerTest', 'broom')) {
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



##################################################################################
# PS5: LMMs and GLMMs
##################################################################################

# Questions (which are repeated below:

# 1. Write out the two covariance matrices given by the model lmer(Y ~ X + (1 | ROUTE)) 
# assuming that there are 3 routes and 3 stations per route. Do the same for the model 
# lmer(Y ~ 1 + (X | ROUTE)). [This second one is harder. This is also not a model you 
# are likely to use, but it is easier to work out than a model with a random effects 
# term for X that you would likely use, like lmer(Y ~ 1 + (1 + X | ROUTE)).]

# 2. Investigate the distribution of the LM, GLM, LMM, and GLMM estimators of b1. Do 
# they estimate the same thing? Are any of the estimators biased? What happens as you 
# increase the route-to-route variation (increase sd.RUGR)? [NOTE: See the pdfs that I 
# already produced.]

# 3. Compare the power curves of LM, GLM, LMM, and GLMM with the same data sets with 
# sd.RUGR = 0.1 and 1.2. Which estimator would you pick? How does the value of 
# sd.RUGR affect type I errors? [NOTE: See the pdfs that I already produced.][NOTE: 
# The black and blue lines sometimes coincide.]

##################################################################################
# Input and organize the grouse data for RUGR and PERIOD 3
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





# -----------------
# 1. Write out the two covariance matrices given by the model lmer(Y ~ X + (1 | ROUTE)) 
# assuming that there are 3 routes and 3 stations per route. Do the same for the model 
# lmer(Y ~ 1 + (X | ROUTE)). [This second one is harder. This is also not a model you 
# are likely to use, but it is easier to work out than a model with a random effects 
# term for X that you would likely use, like lmer(Y ~ 1 + (1 + X | ROUTE)).]
# -----------------

# for lmer(Y ~ X + (1 | ROUTE))
# n by n matrix of covariances, which would be 'blocked' by routes.


# lmer(Y ~ 1 + (X | ROUTE))
# ??



##########################################
# Logitnormal-Binomial for observations at the STATION level


# Simulate data w with a logitnormal-binomial at the station level
lnorm_b_sim <- function(groups, xs, b0, b1, sd) {
    # Inner function to do a single group's simulation
    group_sim <- function(group_i) {
        xs_i <- xs[groups == group_i]
        n <- length(xs_i)
        E <- rnorm(n = 1, mean = 0, sd = sd)
        prob_obs <- inv_logit(b0 + (b1 * xs_i) + E)
        rbinom(n = n, size = 1, prob = prob_obs)
    }
    
    if (class(groups) == 'factor') {
        unq_groups <- levels(groups)
    } else {
        unq_groups <- unique(groups)
    }
    
    c(lapply(unq_groups, group_sim), recursive = TRUE)
}



comp_methods <- function(x_vec, me_vec, sim_fun, ...) {
    
    sim_y <- sim_fun(...)
    
    models <- list(
        lm = lm(sim_y ~ x_vec),
        glm = glm(sim_y ~ x_vec, family = "binomial"),
        lmm = lmer(sim_y ~ x_vec + (1 | me_vec), 
                   control = lmerControl(calc.derivs = FALSE)),
        glmm = glmer(sim_y ~ x_vec + (1 | me_vec), family = "binomial", nAGQ = 0,
                     control = glmerControl(calc.derivs = FALSE))
    )
    
    lapply(names(models), 
           function(m){
               data_frame(
                   method = m,
                   b1.est = ifelse(grepl('mm', m), summary(models[[m]])$coef[2,1], 
                                   models[[m]]$coef[2]),
                   P = ifelse(m == 'lmm', summary(models[[m]])$coef[2,5], 
                              summary(models[[m]])$coef[2,4])
               )
           }
    ) %>% bind_rows
}



# -----------------
# 2. Investigate the distribution of the LM, GLM, LMM, and GLMM estimators of b1. 
# Do they estimate the same thing? Are any of the estimators biased? What happens 
# as you increase the route-to-route variation (increase sd.RUGR)? 
# -----------------

# The GLM and GLMM appear to estimate the same thing, as do LM and LMM.
# GLM/GLMM are pretty biased.

# Both mixed effects models seem to *maybe* have greater variability in b1 estimators 
# when sd increases.




b0_true <- -1.9
b1_true <- 0
sd_true <- 0.01

nsims <- 1000


RNGkind("L'Ecuyer-CMRG")
set.seed(333)
s0 <- Sys.time()
comp_sims <- mclapply(1:nsims, 
                      function(i){
                          comp_methods(g3$WINDSPEEDSQR, g3$ROUTE, lnorm_b_sim, 
                                       g3$ROUTE, g3$WINDSPEEDSQR,
                                       b0_true, b1_true, sd_true)
                      }, 
                      mc.cores = ncpus) %>% 
    bind_rows
s1 <- Sys.time()
s1 - s0
# Time difference of 3.698444 mins





# This code takes a long time, to just see the pdfs in the Problem Sets Folder
# pdf("PS5 Q2 Estimators of b1 with b0=-1.9, b1=0, sd.RUGR=0.01.pdf", width=4, height=8)
# par(mfrow = c(4,1))

# hist(output$b1.est.lm, main = "LM estimates of b1", freq=F, breaks=40)
# lines(.001*(-1000:1000), dnorm(.001*(-1000:1000), mean=mean(output$b1.est.lm), sd=sd(output$b1.est.lm)), col="red")

# hist(output$b1.est.glm, main = "GLM estimates of b1", freq=F, breaks=40)
# lines(.01*(-1000:1000), dnorm(.01*(-1000:1000), mean=mean(output$b1.est.glm), sd=sd(output$b1.est.glm)), col="red")

# hist(output$b1.est.lmm, main = "LMM estimates of b1", freq=F, breaks=40)
# lines(.001*(-1000:1000), dnorm(.001*(-1000:1000), mean=mean(output$b1.est.lmm), sd=sd(output$b1.est.lmm)), col="red")

# hist(output$b1.est.glmm, main = "GLMM estimates of b1", freq=F, breaks=40)
# lines(.01*(-1000:1000), dnorm(.01*(-1000:1000), mean=mean(output$b1.est.glmm), sd=sd(output$b1.est.glmm)), col="red")
# dev.off()




# -----------------
# 3. Compare the power curves of LM, GLM, LMM, and GLMM with the same data sets with 
# sd.RUGR = 0.1 and 1.2. Which estimator would you pick? How does the value of sd.RUGR 
# affect type I errors? [NOTE: The black and blue lines sometimes coincide.]
# -----------------

# In the power tests with lowwer sd.RUGR, all models are nearly exactly the same. The
# GLMM might have slightly elevated type I error, as shown by it's higher value when
# b1 = 0.

# In the power tests with higher sd.RUGR:
# "Power" is highest for the non-mixed models, but I use quotes bc it's not valid power
# when type I error is out of control. At b1 = 0, both LM and GLM show elevated type I
# error. Both LMM and GLMM have much better type I error control, but the GLMM seems a 
# little elevated compared to the LMM. If GLMM's type I error control is actually okay,
# then it seems to be a more powerful type of model.


b0_true <- -1.9
sd_true <- 1.2

b1_range <- seq(0, 0.5, 0.1)

RNGkind("L'Ecuyer-CMRG")
set.seed(111)
s0 <- Sys.time()
power_sims <- mclapply(
    b1_range, 
    function(b1_i){
        replicate(nsims, 
                  comp_methods(g3$WINDSPEEDSQR, g3$ROUTE, lnorm_b_sim, 
                               g3$ROUTE, g3$WINDSPEEDSQR,
                               b0_true, b1_i, sd_true), 
                  simplify = FALSE) %>% 
            bind_rows
        }, mc.cores = ncpus) %>% 
    bind_rows
s1 <- Sys.time()
s1 - s0
# Time difference of 22.83665 mins (on my machine)

# If using > 6 cores
RNGkind("L'Ecuyer-CMRG")
set.seed(111)
s0 <- Sys.time()
power_sims <- mclapply(
    1:nsims, 
    function(i){
        lapply(b1_range, 
               function(b1_i){
                   comp_methods(g3$WINDSPEEDSQR, g3$ROUTE, lnorm_b_sim, 
                                g3$ROUTE, g3$WINDSPEEDSQR,
                                b0_true, b1_i, sd_true)
               }) %>% 
            bind_rows
    }, mc.cores = ncpus) %>% 
    bind_rows
s1 <- Sys.time()
s1 - s0




# Plot power curves
# This code takes a long time, to just see the pdfs in the Problem Sets Folder

# output$rejected.lm <- output$P.lm < 0.05
# power.lm <- aggregate(output$rejected.lm, by = list(output$b1.true), FUN = mean)
# names(power.lm) <- c("b1", "rejected")

# output$rejected.glm <- output$P.glm < 0.05
# power.glm <- aggregate(output$rejected.glm, by = list(output$b1.true), FUN = mean)
# names(power.glm) <- c("b1", "rejected")

# output$rejected.glmm <- output$P.glmm < 0.05
# power.glmm <- aggregate(output$rejected.glmm, by = list(output$b1.true), FUN = mean)
# names(power.glmm) <- c("b1", "rejected")

# output$rejected.lmm <- output$P.lmm < 0.05
# power.lmm <- aggregate(output$rejected.lmm, by = list(output$b1.true), FUN = mean)
# names(power.lmm) <- c("b1", "rejected")

# pdf("PS5 Q3 Power for b1 with b0=-1.9, sd.RUGR=0.01.pdf", width=6, height=6)

# par(mfrow=c(1,1))
# plot(rejected ~ b1, data=power.lm, typ="l", main="LM (black), GLM (blue), LMM (green), GLMM (red)", ylim=c(0, max(power.glm$rejected, power.lm$rejected)))
# lines(rejected ~ b1, data=power.glm, col="blue")
# lines(rejected ~ b1, data=power.lmm, col="green")
# lines(rejected ~ b1, data=power.glmm, col="red")
# lines(c(0,10), c(.05,.05), lty=2)

# dev.off()









# # Moving this file to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp Nell_PS5.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/'")
# )
