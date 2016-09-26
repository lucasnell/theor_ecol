
# Install required packages if they're not already installed, then load
for (f in c('dplyr', 'ggplot2', 'grid', 'parallel')) {
    if (!f %in% rownames(installed.packages())) {
        install.packages(f, dependencies = TRUE)
    }
    library(f, character.only = TRUE)
}; rm(f)


b1range <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
nrange <- c(50, 100, 200, 500, 1000)
# Data frame of all combinations of above vectors
n_b1 <- as.tbl(expand.grid(b1 = b1range, n = nrange))

nsims <- 1000

# CHANGE THIS IF YOUR COMPUTER HAS FEWER/MORE CPUs
ncpus <- 3


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



# ==========================================================================
# ==========================================================================

# Simulate a binary dataset

# ==========================================================================
# ==========================================================================


# ------------
# Base functions
# ------------
inv_logit <- function(x){
    exp(x)/(1 + exp(x))
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




# ------------
# Distribution of the estimator of b1
# ------------

# 1. Investigate the distribution of the GLM estimator of b1. Is it biased? What is 
# the distribution? Normal? How does the distribution change with the true (simulated)
# value of b1? 


set.seed(111)
b1_sims <- lapply(c(b1range, 0.8), 
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
    facet_grid(b1.true ~ ., scales = 'free_y') + 
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.y = element_blank(), axis.ticks.y = element_blank())




# Seems mostly normal, but right skewed, especially at higher b1
# Doesn't seem biased




# ------------
# Type I errors
# ------------

# 2. Check the type I errors. Does type I error control depend on n? How would this 
# affect your interpretation of p-values in your analyses?

# Added n=10,20 for comparison with very low sample size
set.seed(222)
typeI <- lapply(c(10, 20, nrange), 
                function(n_i){
                    replicate(nsims, bin_sim(n = n_i), simplify = FALSE) %>%
                        bind_rows %>%
                        mutate(n = n_i)# %>%
                        # filter(abs(b1.est) < 10)
                }
    ) %>% 
    bind_rows



typeI %>% 
    group_by(n) %>% 
    summarize(typeI = mean(P < 0.05, na.rm = TRUE)) %>% 
    ggplot(aes(n, typeI)) + 
    geom_line() + 
    geom_hline(yintercept = 0.05, linetype = 3)


# Type I does appear to change much after n=50, but at n=20 it was 0.022 and 
# at n=10 it was 0.00!

# I'd obviously not feel comfortable saying that smaller sample sizes eliminate Type I
# error. It's probably more likely that at such small sample sizes, our power is 
# incredibly low, so we're not rejecting *anything*.

# 9-21 --> Sample sizes *should not* influence Type I error. Something is wrong with 
# our statistic in this case.



# ------------
# Consistency of b1
# ------------

# Check the GLM estimator of b1 for consistency in the range `nrange`
# Does the variance in the estimator depend on the true value of b1?

set.seed(333)
consist <- lapply(seq(nrow(n_b1)),
                  function(i){
                      n_i <- n_b1$n[i]
                      b1_i <- n_b1$b1[i]
                      replicate(nsims, bin_sim(n = n_i, b1 = b1_i),
                                simplify = FALSE) %>%
                          bind_rows %>%
                          mutate(n = n_i, b1.true = b1_i) %>%
                          filter(abs(b1.est) < 10)
                  }
    ) %>%
    bind_rows


# RNGkind("L'Ecuyer-CMRG")
# set.seed(333)
# consist <- mclapply(seq(nrow(n_b1)),
#                      function(i){
#                          n_i <- n_b1$n[i]
#                          b1_i <- n_b1$b1[i]
#                          replicate(nsims, bin_sim(n = n_i, b1 = b1_i),
#                                    simplify = FALSE) %>%
#                              bind_rows %>%
#                              mutate(n = n_i, b1.true = b1_i) %>%
#                              filter(abs(b1.est) < 10)
#                      },
#                      mc.cores = ncpus) %>%
#     bind_rows
# 
# # My parallelization results:
# 
# # Using standard lapply...
# # user  system elapsed 
# # 101.563   1.154 103.290
# 
# # Using mclapply with 3 cores...
# # user  system elapsed 
# # 88.993   1.093  48.116



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


set.seed(444)
powers <- lapply(seq(nrow(n_b1)),
                  function(i){
                      n_i <- n_b1$n[i]
                      b1_i <- n_b1$b1[i]
                      replicate(nsims, bin_sim(n = n_i, b1 = b1_i, lm_too = TRUE),
                                simplify = FALSE) %>%
                          bind_rows %>%
                          mutate(n = n_i, b1.true = b1_i) %>%
                          filter(abs(b1.est) < 10)
                  }
    ) %>%
    bind_rows

powers$rejected <- powers$P < 0.05

powers_sum <- powers %>%
    group_by(method, n, b1.true) %>%
    summarize(rejected = mean(rejected)) %>% 
    ungroup
powers_sum


powers_sum %>% 
    ggplot(aes(n, rejected, color = factor(method))) + 
    geom_line() +
    theme_bw() + 
    facet_grid(. ~ b1.true)


powers_sum %>% 
    ggplot(aes(b1.true, rejected, color = factor(method))) + 
    geom_line() +
    theme_bw() + 
    facet_grid(. ~ n)












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
#           "&& cp Nell_PS3.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/'")
# )




