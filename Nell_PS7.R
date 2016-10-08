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

# formerly rg3
d <- grouse_df %>%
    filter(PERIOD == 3, species == 'Ruffed Grouse')

# formerly rg3_route
w <- d %>% 
    select(detected, ROUTE, STATION, WIND_SPEED, WINDSPEEDSQR) %>% 
    group_by(ROUTE) %>% 
    summarize(
        RUGR = sum(detected),
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




##################################################################################
# ML estimation of the beta-binomial distribution
##################################################################################


# Simulate data w with a beta-binomial at the route level (from PS4)
simulate.w.betabinomial <- function (w, b0, b1, theta) {

	w.sim.Y <- array(-1, dim=dim(w)[1])
	counter <- 1
	for(i in levels(w$ROUTE)){
		n <- w$STATIONS[i == w$ROUTE]
		p <- inv_logit(b0 + b1*(w$WINDSPEEDSQR[counter]))
		probRoute <- rbeta(n=1, shape1=p * theta, shape2=(1-p) * theta)
		w.sim.Y[counter] <- rbinom(n=1, size=n, prob=probRoute)
		counter <- counter + 1
	}
	return(w.sim.Y)
}

# Graph Beta-binomial simulations
par(mfrow=c(3,1))

# Beta-binomial distribution
p <- 0.2
theta <- 5
shape1 <- p * theta
shape2 <- (1-p) * theta
sim.bb <- rbeta(n=1000, shape1=shape1, shape2=shape2)
hist(sim.bb, main='Simulation of beta distribution', xlab="Route probability", freq = F)

# Simulated data
b0 <- logit(p)
b1 <- 0
w.sim <- w
w.sim$RUGR <- simulate.w.betabinomial(w = w, b0 = b0, b1 = b1, theta = theta)

hist(w.sim$RUGR/w.sim$STATIONS, main=paste('Beta-binom: mean = ', .001*round(1000*mean(w.sim$RUGR/w.sim$STATIONS)),'   sd=',.001*round(1000*sd(w.sim$RUGR/w.sim$STATIONS)), sep=''), xlab="RUGR/STATION")

# Real data
hist(w$RUGR/w$STATIONS, main=paste('Real data: mean = ', .001*round(1000*mean(w$RUGR/w$STATIONS)),'   sd=',.001*round(1000*sd(w$RUGR/w$STATIONS)), sep=''), xlab="RUGR/STATION")


##################################################################################
# Probability distribution and likelihood functions

# Probability distribution function for a betabinomial distribution from the 
# library "emdbook"
dbetabinom <- function(y, prob, size, theta, log = FALSE) {
	v <- lfactorial(size) - lfactorial(y) - lfactorial(size - y) - 
	    lbeta(theta * (1 - prob), theta * prob) + 
	    lbeta(size - y + theta * (1 - prob), y + theta * prob)
	if (sum((y%%1) != 0) != 0) {
		warning("non-integer x detected; returning zero probability")
		v[n] <- -Inf
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

# Log-likelihood function for the betabinomial given data Y (vector of successes), Size (vector of number of trials), and independent variable X (WINDSPEEDSQR) with b1 = 0 in terms of parameters prob and theta
dbetabinom_LLF0 <- function(parameters, Y, size, X) {
	theta <- parameters['theta']
	b0 <- parameters['b0']
	
	prob <- inv_logit(b0)
	LL <- -sum(dbetabinom(y = Y, size = size, prob = prob, theta = theta, log = TRUE))

	return(LL)
}

##################################################################################
# Fit to real data

z <- optim(fn = dbetabinom_LLF, par = c(theta = 1, b0 = 0, b1 = 0.1), 
           Y = w$RUGR, size = w$STATIONS, 
           X = w$WINDSPEEDSQR, method = "BFGS")
z$par

# Distribution of the estimators of b0, b1, and theta

b0 <- -1.136343
b1 <- -0.300938
theta <- 5.684785

w.sim <- w

nsims <- 1000
output <- data.frame(b0.est=array(0,dim=nsims), b1.est=0, theta.est.glm=0)
for(i in 1:nsims){
	w.sim$RUGR <- simulate.w.betabinomial(w = w, b0 = b0, b1 = b1, theta = theta)
	z.sim <- optim(fn = dbetabinom_LLF, par = c(theta = 1, b0 = 0, b1 = 0.1), Y = w.sim$RUGR, size = w.sim$STATIONS, X = w.sim$WINDSPEEDSQR, method = "BFGS")
	
	output$b0.est[i] <- z.sim$par['b0']
	output$b1.est[i] <- z.sim$par['b1']
	output$theta.est[i] <- z.sim$par['theta']
}

par(mfrow = c(3,1))

hist(output$b0.est, main = paste('true b0 =', .001*round(1000*b0), ', mean b0.est =', .001*round(1000*mean(output$b0.est))), xlab = "ML estimates of b0", freq=F, breaks=100)
lines(.01*(-1000:1000), dnorm(.01*(-1000:1000), mean=mean(output$b0.est), sd=sd(output$b0.est)), col="red")

hist(output$b1.est, main = paste('true b1 =', .001*round(1000*b1), ', mean b1.est =', .001*round(1000*mean(output$b1.est))), xlab = "ML estimates of b1", freq=F, breaks=100)
lines(.01*(-1000:1000), dnorm(.01*(-1000:1000), mean=mean(output$b1.est), sd=sd(output$b1.est)), col="red")

hist(output$theta.est, main = paste('true theta =', .001*round(1000*theta), ', mean theta.est =', .001*round(1000*mean(output$theta.est))), xlab = "ML estimates of theta", freq=F, breaks=100)
lines(.01*(-1000:1000), dnorm(.01*(-1000:1000), mean=mean(output$theta), sd=sd(output$theta)), col="red")

# 1. Is the ML estimator of b0, b1, or theta biased?

# correlation among parameters
par(mfrow = c(3,1))

plot(output$b1.est ~ output$b0.est, main = 'b1 vs. b0')
plot(output$theta.est ~ output$b0.est, main = 'theta vs. b0')
plot(output$theta.est ~ output$b1.est, main = 'theta vs. b1')

# 2. Explain the pattern of correlation among paramaters. Do these make sense to you?

# Bootstrap testing H0: b1 = 0 from the distribution of the estimator of b1

par(mfrow = c(2,1))

hist(output$b1.est, main = paste('Fraction > 0 = ', .0001*round(10000*mean(output$b1.est > 0)), ', P-value = ', 2 * .0001*round(10000*mean(output$b1.est > 0))), xlab = "ML estimates of b1", freq=F, breaks=100)
lines(c(0, 0), c(0,10), col='red')

# Bootstrap testing H0: b1 = 0 under the null hypothesis

z0 <- optim(fn = dbetabinom_LLF0, par = c(theta = 1, b0 = 0), Y = w$RUGR, size = w$STATIONS, X = w$WINDSPEEDSQR, method = "BFGS")
z0$par

b0 <- -1.541177
b1true <- b1
b1 <- 0
theta <- 5.190277

w.sim <- w

nsims <- 1000
output <- data.frame(b0.est=array(0,dim=nsims), b1.est=0, theta.est.glm=0)
for(i in 1:nsims){
	w.sim$RUGR <- simulate.w.betabinomial(w = w, b0 = b0, b1 = b1, theta = theta)
	z.sim <- optim(fn = dbetabinom_LLF, par = c(theta = 1, b0 = 0, b1 = 0.1), Y = w.sim$RUGR, size = w.sim$STATIONS, X = w.sim$WINDSPEEDSQR, method = "BFGS")
	
	output$b0.est[i] <- z.sim$par['b0']
	output$b1.est[i] <- z.sim$par['b1']
	output$theta.est[i] <- z.sim$par['theta']
}

hist(output$b1.est, main = paste('Fraction < b1true = ', .0001*round(10000*mean(output$b1.est < b1true)), ', P-value = ', 2 * .0001*round(10000*mean(output$b1.est < b1true))), xlab = "ML estimates of b1 under H0", freq=F, breaks=100)
lines(c(b1true, b1true), c(0,10), col='red')

# 3. Compare the two bootstrap methods for determining whether b1 differs from zero. Which method is correct?

# Likelihood ratio test (LRT)
z <- optim(fn = dbetabinom_LLF, par = c(theta = 1, b0 = 0, b1 = 0.1), Y = w$RUGR, size = w$STATIONS, X = w$WINDSPEEDSQR, method = "BFGS")
z$par

z0 <- optim(fn = dbetabinom_LLF0, par = c(theta = 1, b0 = 0), Y = w$RUGR, size = w$STATIONS, X = w$WINDSPEEDSQR, method = "BFGS")
z0$par

logLik <- -z$value
logLik0 <- -z0$value

pchisq(2*(logLik - logLik0), df = 1, lower.tail = F)

# 4. Compare the bootstrap tests of H0:b1 = 0 to the LRT. Which do you think gives the best test? What would you want to do to validate the LRT? (Bonus points if you go ahead and do this validation!)

##################################################################################
# Comparison among all methods
##################################################################################

# 5. Perform all of the analyses that we have done so far to test whether WINDSPEEDSQR has an effect on detecting RUGR in PERIOD 3. What are the pros and cons for each analysis? I have given some of the code below that you could use for ML estimation of the beta-binomial and bootstrapping the glmms.

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

library(lme4)
library(lmerTest)
require("boot")

# This is the inverse logit function
inv_logit <- function(x){
	1/(1 + exp(-x))
}

# Simulate data w with a beta-binomial at the route level (from PS4)
simulate.w.betabinomial <- function (w, b0, b1, theta) {

	w.sim.Y <- array(-1, dim=dim(w)[1])
	counter <- 1
	for(i in levels(w$ROUTE)){
		n <- w$STATIONS[i == w$ROUTE]
		p <- inv_logit(b0 + b1*(w$WINDSPEEDSQR[counter]))
		probRoute <- rbeta(n=1, shape1=p * theta, shape2=(1-p) * theta)
		w.sim.Y[counter] <- rbinom(n=1, size=n, prob=probRoute)
		counter <- counter + 1
	}
	return(w.sim.Y)
}

# Probability distribution function for a betabinomial distribution from the library "emdbook"
dbetabinom <- function(y, prob, size, theta, log = FALSE) {
	v <- lfactorial(size) - lfactorial(y) - lfactorial(size - y) - lbeta(theta * (1 - prob), theta * prob) + lbeta(size - y + theta * (1 - prob), y + theta * prob)
	if (sum((y%%1) != 0) != 0) {
		warning("non-integer x detected; returning zero probability")
		v[n] <- -Inf
	}
	if (log == T) {
		return(v)
	}else {
		return(exp(v))
	}
}

# Log-likelihood function for the betabinomial given data Y (vector of successes), Size (vector of number of trials), and independent variable X (WINDSPEEDSQR) in terms of parameters prob and theta
dbetabinom_LLF <- function(parameters, Y, size, X) {
	theta <- parameters['theta']
	b0 <- parameters['b0']
	b1 <- parameters['b1']

	prob <- inv_logit(b0 + b1 * X)	
	LL <- -sum(dbetabinom(y = Y, size = size, prob = prob, theta = theta, log = TRUE))

	return(LL)
}

# Log-likelihood function for the betabinomial given data Y (vector of successes), Size (vector of number of trials), and independent variable X (WINDSPEEDSQR) with b1 = 0 in terms of parameters prob and theta
dbetabinom_LLF0 <- function(parameters, Y, size, X) {
	theta <- parameters['theta']
	b0 <- parameters['b0']
	
	prob <- inv_logit(b0)
	LL <- -sum(dbetabinom(y = Y, size = size, prob = prob, theta = theta, log = TRUE))

	return(LL)
}

# Simulate data w with a logitnormal-binomial at the route level
simulate.w.logitnormalbinomial <- function (w, b0, b1, sd) {

	w.sim.Y <- array(-1, dim=dim(w)[1])
	counter <- 1
	for(i in levels(w$ROUTE)){
		nn <- w$STATIONS[i == w$ROUTE]
		e <- rnorm(n = 1, mean = 0, sd = sd)
		probRoute <- inv_logit(b0 + b1*(w$WINDSPEEDSQR[counter]) + e)
		w.sim.Y[counter] <- rbinom(n=1, size=nn, prob=probRoute)
		counter <- counter + 1
	}
	return(w.sim.Y)
}

# Simulate data d with a logitnormal-binomial at the station level (from PS4)
simulate.d.logitnormalbinomial <- function (d, b0, b1, sd) {

	d.sim.Y <- array(-1, dim=dim(d)[1])
	for(i in levels(d$ROUTE)){
		dd <- d[d$ROUTE == i,]
		nn <- dim(dd)[1]
		probRoute <- rnorm(n=1, mean=0, sd=sd)
		probObs <- inv_logit(b0 + b1*dd$WINDSPEEDSQR + probRoute)
		d.sim.Y[d$ROUTE == i] <- rbinom(n=nn, size=1, prob=probObs)
	}
	return(d.sim.Y)
}

##################################################################################
# ROUTE-level data

summary(w)

##################################################################################
# STATION-level data

summary(d)
