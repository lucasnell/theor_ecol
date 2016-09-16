
# Install required packages if they're not already installed, then load
for (f in c('dplyr', 'ggplot2', 'grid')) {
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


# ------------
# Simulating data
# ------------
cont_sim <- function(b0 = 0, b1 = 0, e.sd = 1, n = 10, single_ex = FALSE) {
    X <- rnorm(n = n)
    Y <- b0 + b1*X + rnorm(n = n, mean = 0, sd = e.sd)
    if (single_ex) {
        return(data_frame(x = X, y = Y))
    }
    z <- lm(Y ~ X, data = NULL)
    out_df <- data_frame(
        b0.est = z$coef[1],
        b1.est = z$coef[2],
        resid2.est = mean(z$resid^2),  # "residual sum of squares"
        sigma.est = mean(z$resid^2)*n/(n-2),  # "estimated variance of the random error"
        P = summary(z)$coef[2,4]
    )
    return(out_df)
}


nsims <- 1000

set.seed(111)
sim_df <- replicate(nsims, cont_sim(), simplify = FALSE) %>% bind_rows
sing_df <- cont_sim(single_ex = TRUE)


lm_p <- ggplot(sing_df, aes(x, y)) +
    ggtitle("Single example") + 
    geom_point() +
    geom_smooth(method = 'lm', se = FALSE)

sim_b1_h <- ggplot(sim_df, aes(b1.est)) + 
    geom_histogram(bins = 20, fill = 'dodgerblue') +
    ggtitle(paste('mean =', sprintf("%.3f", mean(sim_df$b1.est)),
                  ' sd =', sprintf("%.3f", sd(sim_df$b1.est))))

grid.newpage()
grid.draw(rbind(ggplotGrob(lm_p), 
                ggplotGrob(sim_b1_h), 
                size = "last"))



# ------------
# Two estimates of sigma
# ------------

# Residual sum of squares
sim_resid2_h <- ggplot(sim_df, aes(resid2.est)) + 
    geom_histogram(bins = 30, fill = 'red') +
    geom_vline(aes(xintercept = 1^2), linetype= 3) +
    ggtitle(paste('mean =', sprintf("%.3f", mean(sim_df$resid2.est)),
                  ' sd =', sprintf("%.3f", sd(sim_df$resid2.est)))) + 
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 150))

# Sigma: estimated variance of the random error
sim_sigma_h <- ggplot(sim_df, aes(sigma.est)) + 
    geom_histogram(bins = 30, fill = 'blue') +
    geom_vline(aes(xintercept = 1^2), linetype= 3) +
    ggtitle(paste('mean =', sprintf("%.3f", mean(sim_df$sigma.est)),
                  ' sd =', sprintf("%.3f", sd(sim_df$sigma.est)))) + 
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 150))

grid.newpage()
grid.draw(rbind(ggplotGrob(sim_resid2_h), 
                ggplotGrob(sim_sigma_h), 
                size = "last"))

# ------------
# Bias
# ------------
# Using Mean Signed Deviation (MSD)

# In estimators of sigma
sum((sim_df$resid2.est - 1) / nrow(sim_df))
sum((sim_df$sigma.est - 1) / nrow(sim_df))
# sigma.est is *less* biased

# In b1
sum((sim_df$b1.est - 0) / nrow(sim_df))


# ------------
# Efficiency
# ------------
# Relative
{sum((sim_df$resid2.est - 1)^2) / nrow(sim_df)} / 
{sum((sim_df$sigma.est - 1)^2) / nrow(sim_df)}
# sigma.est is *less* efficient


# ------------
# Type I errors
# ------------
mean(sim_df$P < 0.05)



# ------------
# Consistency
# ------------

nrange <- c(10, 20, 50, 100, 200, 500, 1000)

mult_ns <- data_frame()
set.seed(222)
for(n in nrange){
    output <- replicate(nsims, cont_sim(n = n), simplify = FALSE) %>% 
        bind_rows %>%
        mutate(n = as.integer(n))
    mult_ns <- bind_rows(mult_ns, output)
}; rm(n, output)


consistency <- mult_ns %>% 
    group_by(n) %>% 
    summarize(b1.sd = sd(b1.est))

ggplot(consistency, aes(n, b1.sd)) + geom_line()



# ------------
# Power
# ------------

b1range <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

mult_b1 <- data_frame()
set.seed(333)
for(b1 in b1range){
    output <- replicate(nsims, cont_sim(n = 100, b1 = b1), simplify = FALSE) %>% 
        bind_rows %>%
        mutate(b1.true = b1)
    mult_b1 <- bind_rows(mult_b1, output)
}; rm(b1, output)

mult_b1$rejected <- mult_b1$P < 0.05

powers <- mult_b1 %>% 
    group_by(b1.true) %>% 
    summarize(rejected = mean(rejected))

ggplot(powers, aes(b1.true, rejected)) + geom_line()



# 

mult_np <- data_frame()
set.seed(333)
for(n in nrange){
    output <- replicate(nsims, cont_sim(n = n, b1 = 0.2), simplify = FALSE) %>% 
        bind_rows %>%
        mutate(n = n)
    mult_np <- bind_rows(mult_np, output)
}; rm(n, output)

mult_np$rejected <- mult_np$P < 0.05

powers_n <- mult_np %>% 
    group_by(n) %>% 
    summarize(rejected = mean(rejected))

ggplot(powers_n, aes(n, rejected)) + geom_line()



# ==========================================================================
# ==========================================================================

# Simulate a binary dataset

# ==========================================================================
# ==========================================================================

inv_logit <- function(x){
    exp(x)/(1 + exp(x))
}


bin_sim <- function(b0 = 0, b1 = 0, n = 10) {
    X <- rnorm(n = n)
    Y <- rbinom(n = n, size = 1, prob = inv_logit(b0 + b1*X))
    z <- lm(Y ~ X, data = NULL)
    out_df <- data_frame(
        b0.est = z$coef[1],
        b1.est = z$coef[2],
        resid2.est = mean(z$resid^2),  # "residual sum of squares"
        sigma.est = mean(z$resid^2)*n/(n-2),  # "estimated variance of the random error"
        P = summary(z)$coef[2,4]
    )
    return(out_df)
}

nsims <- 1000

set.seed(111)
sim_df <- replicate(nsims, bin_sim(), simplify = FALSE) %>% bind_rows


# ------------
# Two estimates of sigma
# ------------

# Residual sum of squares
sim_resid2_h <- ggplot(sim_df, aes(resid2.est)) + 
    geom_histogram(bins = 30, fill = 'red') +
    # geom_vline(aes(xintercept = 1^2), linetype= 3) +
    ggtitle(paste('mean =', sprintf("%.3f", mean(sim_df$resid2.est)),
                  ' sd =', sprintf("%.3f", sd(sim_df$resid2.est)))) + 
    coord_cartesian(xlim = c(-0.01, 0.325), ylim = c(0, 150))

# Sigma: estimated variance of the random error
sim_sigma_h <- ggplot(sim_df, aes(sigma.est)) + 
    geom_histogram(bins = 30, fill = 'blue') +
    # geom_vline(aes(xintercept = 1^2), linetype= 3) +
    ggtitle(paste('mean =', sprintf("%.3f", mean(sim_df$sigma.est)),
                  ' sd =', sprintf("%.3f", sd(sim_df$sigma.est)))) + 
    coord_cartesian(xlim = c(-0.01, 0.325), ylim = c(0, 150))

grid.newpage()
grid.draw(rbind(ggplotGrob(sim_resid2_h), 
                ggplotGrob(sim_sigma_h), 
                size = "last"))



# ------------
# Bias
# ------------
# Using Mean Signed Deviation (MSD)

# From simple binomial distribution (does this apply here??)
sigma_known <- 0.5 * 0.5

# In estimators of sigma
sum((sim_df$resid2.est - sigma_known) / nrow(sim_df))
sum((sim_df$sigma.est - sigma_known) / nrow(sim_df))

# In b1
sum((sim_df$b1.est - 0) / nrow(sim_df))


# ------------
# Efficiency
# ------------
# Relative
{sum((sim_df$resid2.est - sigma_known)^2) / nrow(sim_df)} /
{sum((sim_df$sigma.est - sigma_known)^2) / nrow(sim_df)}



# ------------
# Type I errors
# ------------
mean(sim_df$P < 0.05, na.rm = TRUE)

# Not sure why there are any NaNs
sim_df %>% filter(is.na(P))


# ------------
# Consistency
# ------------

nrange <- c(10, 20, 50, 100, 200, 500, 1000)

mult_ns <- data_frame()
set.seed(222)
for(n in nrange){
    output <- replicate(nsims, bin_sim(n = n), simplify = FALSE) %>% 
        bind_rows %>%
        mutate(n = as.integer(n))
    mult_ns <- bind_rows(mult_ns, output)
}; rm(n, output)


consistency <- mult_ns %>% 
    group_by(n) %>% 
    summarize(b1.sd = sd(b1.est))
consistency

ggplot(consistency, aes(n, b1.sd)) + geom_line()



# ------------
# Power
# ------------

b1range <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

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
#     paste("mkdir -p ~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/ps2/' &&",
#           "cd", getwd(),
#           "&& cp Nell_PS2.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/ps2/'")
# )

