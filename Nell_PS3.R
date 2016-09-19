
# Install required packages if they're not already installed, then load
for (f in c('dplyr', 'ggplot2', 'grid')) {
    if (!f %in% rownames(installed.packages())) {
        install.packages(f, dependencies = TRUE)
    }
    library(f, character.only = TRUE)
}; rm(f)

b1range <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
nrange <- c(50, 100, 200, 500, 1000)
n_b1 <- expand.grid(b1 = b1range, n = nrange) %>% as.tbl



##################################################################################
# PS3: Investigating GLMs with simulations
##################################################################################
# Problem Set 3: GLMs

# Please stick with your partner from PS2.

# The code in PS3 first gives a repeat of the application of the LM to the binary 
# simulation data. Then there is an example of fitting a data set using a binary 
# (binomial) GLM. Your job is to adapt the code for LMs for use with GLMs so that you
# can fit simulated data with a binary GLM.

# 1. Investigate the distribution of the GLM estimator of b1. Is it biased? What is 
# the distribution? Normal? How does the distribution change with the true (simulated)
# value of b1? 

# Warning: for smaller n and larger b1, sometimes glm() gives crazy answers. You can 
# get rid of them with the code:

# output <- output[abs(output$b1.est) < 10,]

# which removes values that have values less than -10 or greater than 10. Normally it 
# is a bad thing to throw out estimates like this without checking that it doesn't 
# affect your results, but for now don't worry.

# 2. Check the type I errors. Does type I error control depend on n? How would this 
# affect your interpretation of p-values in your analyses?

# 3. Check the GLM estimator of b1 for consistency in the range

# nrange <- c(50, 100, 200, 500, 1000)

# Does the variance in the estimator depend on the true value of b1?

# 4. Compare the power of GLMs and LMs with the same data sets (i.e., in the loop over
# simulated data sets, test the hypothesis H0: b1 = 0 with both LM and GLM). Is the 
# GLM more or less powerful than the LM? Does the difference depend on the sample size 
# n?





##################################################################################
# Do the following analyses using a binomial GLM

# 1. Distribution of the estimator of b1
# 2. Type I errors
# 3. Consistency
# 4. Power
# 5. Do a comparison of power between the GLM and the LM

# As an example of how to use glm(), try this:

# inv.logit <- function(x){
# 	exp(x)/(1 + exp(x))
# }
# 
# n <- 10
# b0 <- 0
# b1 <- 1
# 
# X <- rnorm(n=n)
# Y <- rbinom(n=n, size=1, prob=inv.logit(b0 + b1*X))
# 
# z <- glm(Y ~ X, family="binomial")
# summary(z)
# 
# par(mfrow=c(1,1))
# plot(Y ~ X)
# lines(inv.logit(b0 + b1*X[order(X)]) ~ X[order(X)], col="red")
# lines(inv.logit(z$coef[1] + z$coef[2]*X[order(X)]) ~ X[order(X)], col="blue")




# ==========================================================================
# ==========================================================================

# Simulate a binary dataset

# ==========================================================================
# ==========================================================================

inv_logit <- function(x){
    exp(x)/(1 + exp(x))
}
logit <- function(x){
    log(x/(1 - x))
}



bin_sim <- function(b0 = 0, b1 = 0, n = 10, lm_too = FALSE) {
    X <- rnorm(n = n)
    Y <- rbinom(n = n, size = 1, prob = inv_logit(b0 + b1*X))
    z <- suppressWarnings(glm(Y ~ X, data = NULL, family = binomial))
    out_df <- data_frame(
        b0.est = z$coef[1],
        b1.est = z$coef[2],
        resid2.est = mean(z$resid^2),  # "residual sum of squares"
        sigma.est = mean(z$resid^2)*n/(n-2),  # "estimated variance of the random error"
        P = summary(z)$coef[2,4]
    )
    if (lm_too) {
        z_lm <- suppressWarnings(lm(Y ~ X, data = NULL))
        out_df$method <- 'glm'
        out_df_lm <- data_frame(
            b0.est = z_lm$coef[1],
            b1.est = z_lm$coef[2],
            resid2.est = mean(z_lm$resid^2),
            sigma.est = mean(z_lm$resid^2)*n/(n-2),
            P = summary(z_lm)$coef[2,4],
            method = 'lm'
        )
        return(bind_rows(out_df, out_df_lm))
    } else {
        return(out_df)
    }
}

nsims <- 1000

set.seed(111)
sim_df <- replicate(nsims, bin_sim(), simplify = FALSE) %>% 
    bind_rows %>% filter(abs(b1.est) < 10)


# ------------
# Distribution of the estimator of b1
# ------------


set.seed(111)
b1_sims <- lapply(b1range, 
       function(b){
           replicate(nsims, bin_sim(b1 = b), simplify = FALSE) %>%
               bind_rows %>%
               mutate(b1.true = b) %>%
               filter(abs(b1.est) < 10)
           }
       ) %>% 
    bind_rows

b1_sims %>% 
    ggplot(aes(x = b1.est, fill = factor(b1.true))) +
    geom_histogram(bins = 100) +
    geom_vline(aes(xintercept = b1.true), linetype = 3) +
    facet_grid(b1.true ~ .)

# Seems mostly normal, but right skewed, especially at higher b1
# Doesn't seem biased



# ------------
# Two estimates of sigma
# ------------

# I believe this is the expected value of sigma here
# np(1-p); where E(p) = 0 and n = 10
sigma_exp <- 10 * 0.5^2

# Residual sum of squares
sim_resid2_h <- ggplot(sim_df, aes(resid2.est)) + 
    geom_histogram(bins = 30, fill = 'red') +
    geom_vline(aes(xintercept = sigma_exp), linetype= 3) +
    ggtitle(paste('mean =', sprintf("%.3f", mean(sim_df$resid2.est)),
                  ' sd =', sprintf("%.3f", sd(sim_df$resid2.est)))) + 
    coord_cartesian(xlim = c(0, 22.5), ylim = c(0, 475))

# Sigma: estimated variance of the random error
sim_sigma_h <- ggplot(sim_df, aes(sigma.est)) + 
    geom_histogram(bins = 30, fill = 'blue') +
    geom_vline(aes(xintercept = sigma_exp), linetype= 3) +
    ggtitle(paste('mean =', sprintf("%.3f", mean(sim_df$sigma.est)),
                  ' sd =', sprintf("%.3f", sd(sim_df$sigma.est)))) + 
    coord_cartesian(xlim = c(0, 22.5), ylim = c(0, 475))

grid.newpage()
grid.draw(rbind(ggplotGrob(sim_resid2_h), 
                ggplotGrob(sim_sigma_h), 
                size = "last"))


# So these are overdispersed?? Weird.





# ------------
# Bias
# ------------
# Using Mean Signed Deviation (MSD)

# In estimators of sigma
sum((sim_df$resid2.est - sigma_exp) / nrow(sim_df))
sum((sim_df$sigma.est - sigma_exp) / nrow(sim_df))

# In b1
sum((sim_df$b1.est - 0) / nrow(sim_df))


# ------------
# Efficiency
# ------------
# Relative
{sum((sim_df$resid2.est - sigma_exp)^2) / nrow(sim_df)} /
{sum((sim_df$sigma.est - sigma_exp)^2) / nrow(sim_df)}



# ------------
# Type I errors
# ------------
mean(sim_df$P < 0.05, na.rm = TRUE)

# 2. Check the type I errors. Does type I error control depend on n? How would this 
# affect your interpretation of p-values in your analyses?


set.seed(333)
typeI <- lapply(nrange, 
                function(n_i){
                    replicate(nsims, bin_sim(n = n_i), simplify = FALSE) %>%
                        bind_rows %>%
                        mutate(n = n_i) %>%
                        filter(abs(b1.est) < 10)
                }
    ) %>% 
    bind_rows



typeI %>% 
    group_by(n) %>% 
    summarize(typeI = mean(P < 0.05, na.rm = TRUE)) %>% 
    ggplot(aes(n, typeI)) + geom_line()


# ------------
# Consistency of b1
# ------------

# Check the GLM estimator of b1 for consistency in the range `nrange`
# Does the variance in the estimator depend on the true value of b1?

set.seed(444)
consist <- lapply(seq(nrow(n_b1)), 
                  function(i){
                      n_i <- n_b1$n[i]
                      b1_i <- n_b1$b1[i]
                      replicate(nsims, bin_sim(n = n_i, b1 = b1_i), simplify = FALSE) %>%
                          bind_rows %>%
                          mutate(n = n_i, b1.true = b1_i) %>%
                          filter(abs(b1.est) < 10)
                  }
    ) %>% 
    bind_rows


consist %>% 
    group_by(n, b1.true) %>% 
    summarize(b1.sd = sd(b1.est)) %>% 
    ggplot(aes(n, b1.sd, color = factor(b1.true))) + geom_line()





# ------------
# Power
# ------------


# 4. Compare the power of GLMs and LMs with the same data sets (i.e., in the loop over
# simulated data sets, test the hypothesis H0: b1 = 0 with both LM and GLM). Is the 
# GLM more or less powerful than the LM? Does the difference depend on the sample size 
# n?



mult_b1 <- data_frame()
set.seed(333)
for(b1 in b1range){
    output <- replicate(nsims, bin_sim(n = 100, b1 = b1), simplify = FALSE) %>% 
        bind_rows %>%
        mutate(b1.true = b1)
    mult_b1 <- bind_rows(mult_b1, output)
}; rm(b1, output)

mult_b1$rejected <- mult_b1$P < 0.05

powers <- mult_b1 %>% 
    group_by(b1.true) %>% 
    summarize(rejected = mean(rejected))
powers

ggplot(powers, aes(b1.true, rejected)) + geom_line()





