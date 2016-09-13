
# Install required packages if they're not already installed, then load
for (f in c('dplyr', 'ggplot2', 'tidyr', 'broom')) {
    if (!f %in% rownames(installed.packages())) {
        install.packages(f, dependencies = TRUE)
    }
    library(f, character.only = TRUE)
}; rm(f)




# =========================
# Assignment:
# =========================
# 1. Use the code in PS2 to examine the inference properties of LMs estimated using 
# Least Squares when applied to simulated regression data with continuous values of Y. 
# Investigate the properties of bias, efficiency, consistency, type I errors, and 
# power for the estimator of the slope b1. Also investigate bias and efficiency of 
# the estimator of the variance of the residuals; is the unbiased estimator the more 
# efficient?

# 2. Use the code in PS2 to examine the inference properties of LMs estimated using 
# Least Squares when applied to simulated regression data with binary values of Y. 
# Investigate the properties of bias, efficiency, consistency, type I errors, and 
# power for the estimator of the slope b1. Also investigate bias and efficiency of 
# the estimator of the variance of the residuals.

# 3. (optional) Derive mathematically the formula for the estimator of LS regression 
# coefficients.





# ========================================================================================
# ========================================================================================

# PS2: Investigating LMs with simulations

# ========================================================================================
# ========================================================================================

b0 <- 0
b1 <- 0
e.sd <- 1
n <- 10

nsims <- 1000
output <- data.frame(b0.est=array(0,dim=nsims), b1.est=0, resid2.est=0, sigma.est=0, P=0)
for(i in 1:nsims){
    X <- rnorm(n=n)
    Y <- b0 + b1*X + rnorm(n=n, mean=0, sd=e.sd)
    z <- lm(Y ~ X)
    output$b0.est[i] <- z$coef[1]
    output$b1.est[i] <- z$coef[2]
    output$resid2.est[i] <- mean(z$resid^2)
    output$sigma.est[i] <- mean(z$resid^2)*n/(n-2)
    output$P[i] <- summary(z)$coef[2,4]
}
par(mfrow = c(2,1))

plot(Y ~ X, main="Single Example")
lines(z$fitted.values[order(X)] ~ X[order(X)])

hist(output$b1.est, main = paste('mean = ', .001*round(1000*mean(output$b1.est)),'   sd=',.001*round(1000*sd(output$b1.est)), sep=''))

# Two estimates of sigma
par(mfrow = c(2,1))

hist(output$resid2, main = paste('mean = ', .001*round(1000*mean(output$resid2)),'   sd=',.001*round(1000*sd(output$resid2)), sep=''))

hist(output$sigma.est, main = paste('mean = ', .001*round(1000*mean(output$sigma.est)),'   sd=',.001*round(1000*sd(output$sigma.est)), sep=''))

# Type I errors
mean(output$P < 0.05)

# Consistency
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
    Y <- b0 + b1*X + rnorm(n=n, mean=0, sd=e.sd)
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

# Power
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
    Y <- b0 + b1*X + rnorm(n=n, mean=0, sd=e.sd)
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


















# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================



# # Moving this file to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp Nell_PS2.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/ps1/'")
# )

