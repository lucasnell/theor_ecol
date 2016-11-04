
# Install required packages if they're not already installed, then load
for (f in c('magrittr', 'dplyr', 'readr', 'tidyr', 'ggplot2', 'parallel', 'lme4')) {
    if (!f %in% rownames(installed.packages())) { install.packages(f) }
    library(f, character.only = TRUE)
}; rm(f)


# PS_Minoow_24Oct16.R

# The first part of this problem set investigates AR(1) and AR(2) proecesses. The 
# second part gives the code that I presented in class that does initial data 
# processing for minnow data. There are


#########################################################################################
# AR(1)

a1 <- -1.05


t.max <- 100
X <- array(0, dim = t.max)
for(t in 2:t.max){
	X[t] <- a1 * X[t-1] + rnorm(1)
}

layout(mat=matrix(c(1,2,1,3), nrow = 2, ncol = 2))
plot(X, typ='l', xlab='Time', main=paste('AR(1) with a1 =', a1, '   lambda =', 
                                         a1, '   variance =', 
                                         format(var(X), digits = 3)))
plot(X[2:t.max] ~ X[1:(t.max-1)], typ='p', xlab='x(t-1)', ylab='x(t)')
lines(c(-20,20), c(-20,20), col='red')
z <- lm(X[2:t.max] ~ X[1:(t.max-1)])
lines(c(-20,20), z$coef[1] + z$coef[2] * c(-20,20), col='blue')
pacf(X, 5, main='')

# Question 1: (a) How does a1 affect the partial autocorrelation? (b) How does a1 
# affect the variation of the stationary distribution? (c) What happens when a1 > 1? 
# when a1 < -1?

# (a) The partial autocorrelation of lag 1 largely mirrors it.

# (b) Higher a1 increases variance.

# (c) For a1 > 1, x increases exponentially from the mean.
#     For a1 < -1, x has increasingly large peaks around the mean (peak amplitudes 
#       increase exponentially).


#########################################################################################
# AR(2)

a1 <- -0.25
a2 <- -0.5

t.max <- 100
X <- array(0, dim = t.max)
for(t in 3:t.max){
	X[t] <- a1 * X[t-1] + a2 * X[t-2] + rnorm(1)
}

B <- matrix(c(a1, a2,1,0), byrow = T, ncol = 2)
lambda <- eigen(B)$values[1]

layout(mat=matrix(c(1,2,1,3), nrow=2, ncol=2))
plot(X, typ='l', xlab='Time', main=paste('AR(1) with a1 =', a1, ' a2 =', a2, 
                                         '   lambda =', format(lambda, digits=2), 
                                         '   ||lambda|| =', 
                                         format(abs(lambda), digits=2), '   variance =',
                                         format(var(X), digits = 3)))
plot(X[2:t.max] ~ X[1:(t.max-1)], typ='l', xlab='x(t-1)', ylab='x(t)')
pacf(X, 5, main='')

# Question 2: (a) How do a1 and a2 affect the partial autocorrelation? (b) How does 
# ||lambda|| affect the variation of the stationary distribution? (c) What happens to 
# the pacf when lambda is a complex number? (d) What happens with ||lambda|| > 1? 

# (a) From 0-0.5, the lag(1) and lag(2) mirror the a1 and a2 approximately. 
#     After sum(c(a1, a2)) > 1, all the autocorrelation seems to be in the lag(1)

# (b) When lambda is a complex number, the relationship between x(t) and x(t-1) becomes
#     cyclical.

# (c) PACF is negative??

# (d) When lambda is > 1, X move exponentially away from the mean



# AR(1) is either boom and bust peaks, or sudden change
# AR(2) can be cyclical



#########################################################################################
# Minnow data
#########################################################################################

options(stringsAsFactors=F)
d <- read.csv(file="./Pine_dat/catch_data.csv", header=T) %>% 
    # filter(isolated_pools == 0) %>% 
    select(-isolated_pools)
f <- read.csv(file="./Pine_dat/flow_data.csv", header=T)
# summary(d)
# summary(f)

# -------------
# Column descriptions
# -------------
sapply(d, class)
#    year       reach     station    rm_start       count        area
# "integer" "character" "character"   "numeric"   "integer"   "numeric"

# year: year sampled
# reach: region name that distinguishes between river areas divided by dams
# station: individual station sampled at (this column is not very useful for analyses;
#          I've kept it here mostly for my own organizational purposes)
# rm_start: "river mile" where sampling took place
# count: number of minnows caught in seine
# area: number of square meters seined; this will be NA when the site was dry on the
#       date it was surveyed


sapply(f, class)
#     year over_1000 over_2000 over_3000 over_4000 over_5000 over_6000 over_7000
# "integer" "integer" "integer" "integer" "integer" "integer" "integer" "integer"

# year: year flow data is from
# over_X: number of days the river was above X cubic feet per second

###########################################################################

# identify stations
d$station_string <- d$station
d$station <- as.factor(d$rm_start)

# remove dry stations
d <- d[!is.na(d$area) & d$area > 0,]

# CPU
d$CPU <- d$count/d$area

# number of years sampled at each station
plot(d$station)

# make annual data set
x <- data.frame(year = unique(d$year), meanCPU = 0, nstation = 0, meanarea = 0)
count <- 0
for(i in unique(d$year)){
	count <- count + 1
	dd <- d[d$year == i,]
	x$year[count] <- i
	x$meanCPU[count] <- mean(dd$CPU)
	x$nstation[count] <- dim(dd)[1]
	x$meanarea[count] <- mean(dd$area)
}

# get total flows per year
f$totalflow <- rowSums(f[,2:8])

# merge with minnow data
x$totalflow <- f$totalflow[match(x$year, f$year)]

# plot CPU and flow
plot(meanCPU/max(meanCPU) ~ year, data = x, typ="b", pch=1)
lines(totalflow/max(totalflow) ~ year, data = x, col = 'red')

# plot CPU vs. flow
plot(meanCPU ~ totalflow, data = x, typ="b", pch=1)

# do a simple regression
z <- lm(meanCPU ~ totalflow, data = x)

# plot residuals
par(mfrow = c(2,1))
plot(z$resid ~ x$year, typ="b")
hist(z$resid, main="")

# try a log transform
minCPU <- min(x$meanCPU[x$meanCPU > 0])
x$logCPU <- log(x$meanCPU + minCPU)

# plot logCPU vs. flow
plot(meanCPU ~ totalflow, data = x, typ="b", pch=1)
plot(logCPU ~ totalflow, data = x, typ="b", pch=1)

#########################################################################################
# look at residual structure of the data

# do a simple regression
z.log <- lm(logCPU ~ totalflow, data = x)

# plot residuals
par(mfrow = c(2,1))
plot(z.log$resid ~ x$year, typ="b")
hist(z.log$resid, main="")

# do a z-transform for logCPU and U
x$logCPU <- (x$logCPU - mean(x$logCPU))/sd(x$logCPU)
x$U <- x$totalflow
x$U <- (x$U - mean(x$U))/sd(x$U)

# plot logCPU and U
par(mfrow = c(1,1))
plot(logCPU ~ year, data = x, typ = "b", pch = 1, ylim = c(-2,2))
lines(U ~ year, data = x, col = 'red')

# lm
z <- lm(logCPU ~ U, data = x)
summary(z)

par(mfrow = c(2,2))
plot(z$resid ~ x$year, typ="b")
hist(z$resid, main="")
acf(z$resid, main="ACF", lag.max = 6)
pacf(z$resid, main="PACF", lag.max = 6)

#########################################################################################
# Simulate the minnow data
#########################################################################################

# Question 3: Suppose the minnow data are given by the model:

U <- x$U
year <- x$year

a1 <- 0.5
b1 <- 0.649 #from the lm on the real data

t.max <- 22
X <- array(0, dim = t.max)
X[1] <- x$logCPU[1]

for(t in 2:t.max){
	X[t] <- a1 * X[t-1] + b1 * U[t] + rnorm(1)
}

# plot X and U
par(mfrow = c(1,1))
plot(X/max(X) ~ year, typ="b", pch=1)
lines(U/max(U) ~ year, col = 'red')

z.sim <- lm(X ~ U)
summary(z.sim)
z <- lm(logCPU ~ U, data = x)

# plot residuals of X (and x$logCPU in blue), along with histogram, ACF and PACF
par(mfrow = c(2,2))
plot(z.sim$resid ~ year, typ="b")
lines(z$resid ~ year, typ="b", col="red")
hist(z.sim$resid, main="")
acf(z.sim$resid, main="ACF", lag.max = 6)
pacf(z.sim$resid, main="PACF", lag.max = 6)

# How to you think the value of a1 affects the type I error control for the test 
# H0: b1 = 0?


# It appears to increase type I error problems as a1 increases, but it wasn't as 
# dramatic as I anticipated.




U <- x$U
year <- x$year

# a1 <- 0.5
# b1 <- 0.649 #from the lm on the real data
rm(a1, b1)

minnow_sim <- function(a1, b1, nsim = 1000, seed = NULL) {
    one_sim <- function(a1, b1, ...){
        X <- numeric(length(U))
        X[1] <- x$logCPU[1]
        # Have to use for-loop bc epsilon at time t-1 affects the next iteration's X[t-1]
        for(t in 2:length(X)){
            X[t] <- a1 * X[t-1] + b1 * U[t] + rnorm(1)
        }
        z.sim <- lm(X ~ U)
        return(summary(z.sim)$coef[2, 'Pr(>|t|)'])
    }
    if (!is.null(seed)) {set.seed(seed)}
    sapply(1:nsim, one_sim, a1 = a1, b1 = b1)
}


typeI <- data_frame(a1 = rep(seq(0, 1, 0.2), each = 1000))
set.seed(999)
typeI$p <- do.call(c, lapply(seq(0, 1, 0.2), minnow_sim, b1 = 0))

typeI %>% 
    group_by(a1) %>% 
    summarize(rejected = mean(p < 0.05))




# arima(X, order = c(1,0,0), xreg = U)















# 
# # Moving this file to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp Nell_PS_Minnow.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/'")
# )
