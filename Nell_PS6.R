
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
# PS6: LMMs and GLMMs
##################################################################################
# Questions (which are repeated below):

# 1. Compare three methods for analyzing a linear model with 10 routes and 20 stations per route, with routes having the same intercept but different slopes: the methods are (i) 10 separate linear models, (ii) a single linear model with categorical terms for route and route:slope interactions, and (iii) a linear mixed model with route as a random effect for both slope and intercept. For the simulations, set b0.mean <- 1, b0.sd <- 0, b1.mean <- .2, b1.sd <- .05, and e.sd <- 0.1. What are the mean and variance in slopes calculated with each of the three methods?

# 2. Compare two methods for analyzing a linear model with 10 routes and 20 stations per route, with routes having the same slope but different intercepts: the methods are (i) a single linear model with categorical terms for route, and (ii) a linear mixed model with route as a random effect. For the simulations, set b0.mean <- 1, b0.sd <- 1, b1.mean <- .2, b1.sd <- 0, and e.sd <- 1. Note that you will also have to change the models for fitting the data. Which method gives on average the lowest P-values for the significance of the slope? You will have to run the simulations several times for this.

# 3. Compare three methods for analyzing a binary model with 10 routes and 20 stations per route, with routes having different slopes: the methods are (i) 10 separate linear models, (ii) a single linear model with categorical terms for route and route:slope interactions, and (iii) a linear mixed model with route as a random effect for both slope and intercept. For the simulations, set b0.mean <- 0, b0.sd <- 0, b1.mean <- 1, and b1.sd <- .5. What are the mean and variance in slopes calculated with each of the three methods? Are the patterns you get different from those in questions #1 for linear models?

# 4. Compare two methods for analyzing a binary model with 10 routes and 20 stations per route, with routes having the same slope but different intercepts: the methods are (i) a single linear model with categorical terms for route, and (ii) a linear mixed model with route as a random effect. For the simulations, set b0.mean <- 0, b0.sd <- .5, b1.mean <- .2, and b1.sd <- 0. Note that you will also have to change the models for fitting the data. Which method gives on average the lowest P-values for the significance of the slope? You will have to run the simulations several times for this.

# 5. Bonus (optional) question: Write code and check the type I error rate for the glmm testing the hypothesis that there is no variance in the slopes among routes when there is variation in the intercept.

library(lme4)
library(lmerTest)

# Anova() in {car} works better than anova() in {base}
install.packages("car")
library(car)

####################################################################
# NOTE: lmer() likes it better with nstations=20, so I changed it to this. But there are still sometimes error messages. Just ignore them.

# LMM with variation in the intercept and slope
nroutes <- 10
nstations <- 20

b0.mean <- 1
b0.sd <- 0
b1.mean <- .2
b1.sd <- .05
e.sd <- .1

x.mean <- 0
x.sd <- 1

# Simulation of linear model
d <- data.frame(route=array(1:nroutes, dim=nroutes*nstations), X=0, Y=0)
d <- d[order(d$route),]

for(i in 1:nroutes) {
	b0 <- b0.mean + rnorm(n=1, mean=0, sd=b0.sd)
	b1 <- b1.mean + rnorm(n=1, mean=0, sd=b1.sd)
	x <- rnorm(n=nstations, mean=x.mean, sd=x.sd)
	d$X[d$route == i] <- x
	d$Y[d$route == i] <- b0 + b1 * x + rnorm(n=10, mean=0, sd=e.sd)
}

# Plot the data along with separate linear fits for each route
par(mfrow=c(2,1))
plot(Y ~ X, data = d, main="Fits from LMs")
b0.est.list <- 0
b1.est.list <- 0
for(i in 1:nroutes) {
	x <- d$X[d$route == i]
	points(Y[d$route == i] ~ x, data = d, col=i)
	z <- lm(Y[d$route == i] ~ x, data = d)
	lines(z$fitted[order(x)] ~ x[order(x)], col=i)
	b0.est.list[i] <- z$coef[1]
	b1.est.list[i] <- z$coef[2]
}

# Calculate the mean and sd of estimates from the separate LMs
b0.mean.sd <- c(mean(b0.est.list), sd(b0.est.list))
b1.mean.sd <- c(mean(b1.est.list), sd(b1.est.list))
b0.mean.sd
b1.mean.sd

# Estimate LM with a categorical variable for route affecting the slope and the intercept
z.lm <- lm(Y ~ 0 + X + factor(route) + X:factor(route), data = d)
summary(z.lm)
Anova(z.lm)

b0.mean.sd <- c(mean(z.lm$coef[2:11]), sd(z.lm$coef[2:11]))
b1.mean.sd <- c(mean(z.lm$coef[1] + c(0, z.lm$coef[12:20])), sd(c(0, z.lm$coef[12:20])))
b0.mean.sd
b1.mean.sd

# Estimate LMM with route as a random effect for the slope and the intercept
z.lmm <- lmer(Y ~ X + (1 + X | route), data = d)
summary(z.lmm)

plot(Y ~ X, data = d, main="Fits from LMM")
for(i in 1:nroutes) {
	x <- d$X[d$route == i]
	points(Y[d$route == i] ~ x, data = d, col=i)
	xx <- sort(x)
	y <- coef(z.lmm)$route[i,1] + coef(z.lmm)$route[i,2] * xx
	lines(y ~ xx, col=i)
}

# Test the significance of the slope random effect
z0.lmm <- lmer(Y ~ X + (1 | route), data = d)
P.value <- pchisq(2*(logLik(z.lmm) - logLik(z0.lmm)), df=2, lower.tail = F)
P.value

# This does the same in lmerTest
rand(z.lmm)

# 1. Compare three methods for analyzing a linear model with 10 routes and 10 stations per route, with routes having the same intercept but different slopes: the methods are (i) 10 separate linear models, (ii) a single linear model with categorical terms for route and route:slope interactions, and (iii) a linear mixed model with route as a random effect for both slope and intercept. For the simulations, set b0.mean <- 1, b0.sd <- 0, b1.mean <- .2, b1.sd <- .05, and e.sd <- 0.1. What are the mean and variance in slopes calculated with each of the three methods?

# 2. Compare two methods for analyzing a linear model with 10 routes and 10 stations per route, with routes having the same slope but different intercepts: the methods are (i) a single linear model with categorical terms for route, and (ii) a linear mixed model with route as a random effect. For the simulations, set b0.mean <- 1, b0.sd <- 1, b1.mean <- .2, b1.sd <- 0, and e.sd <- 1. Note that you will also have to change the models for fitting the data. Which method gives on average the lowest P-values for the significance of the slope? You will have to run the simulations several times for this.


####################################################################
# Binary GLMM with variation in the intercept and slope

# This is the inverse logit function
inv.logit <- function(x){
	exp(x)/(1 + exp(x))
}

nroutes <- 10
nstations <- 20

b0.mean <- 0
b0.sd <- 0
b1.mean <- 1
b1.sd <- .5

x.mean <- 0
x.sd <- 1

# Simulation of binary data
d <- data.frame(route=array(1:nroutes, dim=nroutes*nstations), X=0, Y=0)
d <- d[order(d$route),]

for(i in 1:nroutes) {
	b0 <- b0.mean + rnorm(n=1, mean=0, sd=b0.sd)
	b1 <- b1.mean + rnorm(n=1, mean=0, sd=b1.sd)
	x <- rnorm(n=nstations, mean=x.mean, sd=x.sd)
	d$X[d$route == i] <- x
	p <- inv.logit(b0 + b1 * x)
	d$Y[d$route == i] <- rbinom(n = nstations, size = 1, prob = p)
}

# Plot the data along with separate generalized linear fits for each route
par(mfrow=c(2,1))
plot(Y ~ X, data = d, main="Fits from GLMs")
b0.est.list <- 0
b1.est.list <- 0
for(i in 1:nroutes) {
	x <- d$X[d$route == i]
	points(Y[d$route == i] ~ x, data = d, col=i)
	z <- glm(Y[d$route == i] ~ x, family = "binomial", data = d)
	
	xx <- .1*((0:200) - 100)
	y <- inv.logit(z$coef[1] + z$coef[2] * xx)
	lines(y ~ xx, col=i)
	
	b0.est.list[i] <- z$coef[1]
	b1.est.list[i] <- z$coef[2]
}

# Calculate the mean and sd of estimates from the separate GLMs
b0.mean.sd <- c(mean(b0.est.list), sd(b0.est.list))
b1.mean.sd <- c(mean(b1.est.list), sd(b1.est.list))
b0.mean.sd
b1.mean.sd

# Estimate GLM with a categorical variable for route affecting the slope and the intercept
z.glm <- glm(Y ~ 0 + X + factor(route) + X:factor(route), family = "binomial", data = d)
summary(z.glm)
Anova(z.glm)

b0.mean.sd <- c(mean(z.glm$coef[2:11]), sd(z.glm$coef[2:11]))
b1.mean.sd <- c(mean(z.glm$coef[1] + c(0, z.glm$coef[12:20])), sd(c(0, z.lm$coef[12:20])))
b0.mean.sd
b1.mean.sd

# Estimate GLMM with route as a random effect for the slope and the intercept
z.glmm <- glmer(Y ~ X + (1 + X | route), family = "binomial", data = d, glmerControl(calc.derivs=F))
summary(z.glmm)

plot(Y ~ X, data = d, main="Fits from GLMM")
for(i in 1:nroutes) {
	x <- d$X[d$route == i]
	points(Y[d$route == i] ~ x, data = d, col=i)
	xx <- .1*((0:200) - 100)
	y <- inv.logit(coef(z.glmm)$route[i,1] + coef(z.glmm)$route[i,2] * xx)
	lines(y ~ xx, col=i)
}

# Test the significance of the slope random effect
z0.glmm <- glmer(Y ~ X + (1 | route), family = "binomial", data = d, glmerControl(calc.derivs=F))
P.value <- pchisq(2*(logLik(z.glmm) - logLik(z0.glmm)), df=2, lower.tail = F)
P.value

# 3. Compare three methods for analyzing a binary model with 10 routes and 20 stations per route, with routes having different slopes: the methods are (i) 10 separate linear models, (ii) a single linear model with categorical terms for route and route:slope interactions, and (iii) a linear mixed model with route as a random effect for both slope and intercept. For the simulations, set b0.mean <- 0, b0.sd <- 0, b1.mean <- 1, and b1.sd <- .5. What are the mean and variance in slopes calculated with each of the three methods? Are the patterns you get different from those in questions #1 for linear models?

# 4. Compare two methods for analyzing a binary model with 10 routes and 20 stations per route, with routes having the same slope but different intercepts: the methods are (i) a single linear model with categorical terms for route, and (ii) a linear mixed model with route as a random effect. For the simulations, set b0.mean <- 0, b0.sd <- .5, b1.mean <- .2, and b1.sd <- 0. Note that you will also have to change the models for fitting the data. Which method gives on average the lowest P-values for the significance of the slope? You will have to run the simulations several times for this.

# 5. Bonus (optional) question: Write code and check the type I error rate for the glmm testing the hypothesis that there is no variance in the slopes among routes when there is variation in the intercept.
