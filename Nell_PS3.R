
# Install required packages if they're not already installed, then load
for (f in c('dplyr', 'ggplot2', 'grid')) {
    if (!f %in% rownames(installed.packages())) {
        install.packages(f, dependencies = TRUE)
    }
    library(f, character.only = TRUE)
}; rm(f)





##################################################################################
# PS3: Investigating GLMs with simulations
##################################################################################
# Problem Set 3: GLMs

# Please stick with your partner from PS2.

# The code in PS3 first gives a repeat of the application of the LM to the binary simulation data. Then there is an example of fitting a data set using a binary (binomial) GLM. Your job is to adapt the code for LMs for use with GLMs so that you can fit simulated data with a binary GLM.

# 1. Investigate the distribution of the GLM estimator of b1. Is it biased? What is the distribution? Normal? How does the distribution change with the true (simulated) value of b1? 

# Warning: for smaller n and larger b1, sometimes glm() gives crazy answers. You can get rid of them with the code:

# output <- output[abs(output$b1.est) < 10,]

# which removes values that have values less than -10 or greater than 10. Normally it is a bad thing to throw out estimates like this without checking that it doesn't affect your results, but for now don't worry.

# 2. Check the type I errors. Does type I error control depend on n? How would this affect your interpretation of p-values in your analyses?

# 3. Check the GLM estimator of b1 for consistency in the range

# nrange <- c(50, 100, 200, 500, 1000)

# Does the variance in the estimator depend on the true value of b1?

# 4. Compare the power of GLMs and LMs with the same data sets (i.e., in the loop over simulated data sets, test the hypothesis H0: b1 = 0 with both LM and GLM). Is the GLM more or less powerful than the LM? Does the difference depend on the sample size n?



##################################################################################
# Below is first the analysis of binary data using LMs. I will go through this in class on Monday, 19 Sep

inv.logit <- function(x){
	exp(x)/(1 + exp(x))
}

n <- 10
b0 <- 0
b1 <- 1

X <- rnorm(n=n)
Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*X))

par(mfrow=c(1,2))
plot(Y ~ X)
lines(inv.logit(b0 + b1*X[order(X)]) ~ X[order(X)], col="red")
boxplot(X ~ Y, horizontal=TRUE)

# Analyses
n <- 30
b0 <- 0
b1 <- 2

nsims <- 1000
output <- data.frame(b0.est=array(0,dim=nsims), b1.est=0, resid2.est=0, sigma.est=0, P=0, varX=0)
for(i in 1:nsims){
	X <- rnorm(n=n)
	Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*X))
	z <- lm(Y ~ X)
	output$b0.est[i] <- z$coef[1]
	output$b1.est[i] <- z$coef[2]
	output$resid2.est[i] <- mean(z$resid^2)
	output$sigma.est[i] <- mean(z$resid^2)*n/(n-2)
	output$P[i] <- summary(z)$coef[2,4]
	output$varX[i] <- var(X)
}

# 1. Distribution of the estimator of b1
par(mfrow = c(2,1))

plot(Y ~ X, main="Single Example")
lines(z$fitted.values[order(X)] ~ X[order(X)])
lines(inv.logit(b0 + b1*X[order(X)]) ~ X[order(X)], col="red")

hist(output$b1.est, main = paste('mean = ', .001*round(1000*mean(output$b1.est)),'   sd=',.001*round(1000*sd(output$b1.est)), sep=''), freq=F, breaks=40)
lines(.025*(-40:40), dnorm(.025*(-40:40), mean=mean(output$b1.est), sd=sd(output$b1.est)), col="red")

# 2. Type I errors
mean(output$P < 0.05)

# 3. Consistency
b0 <- 0
b1 <- 0
e.sd <- 1
nrange <- c(10, 20, 50, 100, 200, 500, 1000)

nsims <- 1000
output <- data.frame(n=array(0,dim=nsims*length(nrange)), b0.est=0, b1.est=0, resid2.est=0, sigma.est=0, P=0)

counter <- 0
for(n in nrange) for(i in 1:nsims){
	counter <- counter + 1
	X <- rnorm(n=n)
	Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*X))
	z <- lm(Y ~ X)
	output$n[counter] <- n
	output$b0.est[counter] <- z$coef[1]
	output$b1.est[counter] <- z$coef[2]
	output$resid2.est[counter] <- mean(z$resid^2)
	output$sigma.est[counter] <- mean(z$resid^2)*n/(n-2)
	output$P[counter] <- summary(z)$coef[2,4]
}
consistency <- aggregate(output$b1.est, by = list(output$n), FUN = sd)
names(consistency) <- c("n", "b1.sd")

par(mfrow=c(1,1))
plot(b1.sd ~ n, data=consistency, typ="l")

# 4. Power
b0 <- 0
e.sd <- 1
b1range <- c(0, .1, .2, .3, .4, .5)
n <- 100

nsims <- 1000
output <- data.frame(b1.true=array(0,dim=nsims*length(b1range)), b0.est=0, b1.est=0, resid2.est=0, sigma.est=0, P=0)

counter <- 0
for(b1 in b1range) for(i in 1:nsims){
	counter <- counter + 1
	X <- rnorm(n=n)
	Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*X))
	z <- lm(Y ~ X)
	output$b1.true[counter] <- b1
	output$b0.est[counter] <- z$coef[1]
	output$b1.est[counter] <- z$coef[2]
	output$resid2.est[counter] <- mean(z$resid^2)
	output$sigma.est[counter] <- mean(z$resid^2)*n/(n-2)
	output$P[counter] <- summary(z)$coef[2,4]
}
output$rejected <- output$P < 0.05
power <- aggregate(output$rejected, by = list(output$b1.true), FUN = mean)
names(power) <- c("b1", "rejected")

par(mfrow=c(1,1))
plot(rejected ~ b1, data=power, typ="l")

##################################################################################
# Do the following analyses using a binomial GLM

# 1. Distribution of the estimator of b1
# 2. Type I errors
# 3. Consistency
# 4. Power
# 5. Do a comparison of power between the GLM and the LM

# As an example of how to use glm(), try this:

inv.logit <- function(x){
	exp(x)/(1 + exp(x))
}

n <- 10
b0 <- 0
b1 <- 1

X <- rnorm(n=n)
Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*X))

z <- glm(Y ~ X, family="binomial")
summary(z)

par(mfrow=c(1,1))
plot(Y ~ X)
lines(inv.logit(b0 + b1*X[order(X)]) ~ X[order(X)], col="red")
lines(inv.logit(z$coef[1] + z$coef[2]*X[order(X)]) ~ X[order(X)], col="blue")

