library(MASS)
library(tidyverse)
library(parallel)
library(grid)
library(car)



# Custom ggplot2 theme
theme_lan <- function(base_size = 10, base_family = 'Helvetica') {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'bold', size = 11),
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid = element_blank()
        )
}


# Inverse logit function
inv_logit <- function(x){
    exp(x)/(1 + exp(x))
}
# This is the logit function
logit <- function(x){
    log(x/(1 - x))
}

# Discrete-time logistic growth for time t+1
log_growth <- function(N_t, r, K){
    new_N <- N_t + r * N_t * (1 - N_t / K)
    # To avoid negative abundances
    new_N <- ifelse(new_N < 0, 0, new_N)
    # Also rounding bc a partial individual makes no sense
    return(matrix(as.integer(new_N), nrow = 1))
}




# Sine curve with limits and set number of phases
# First function is for one x value, and is used inside the 2nd one and for plotting
u_curve_1x <- function(x, time_len, phases = 1, c_max = 1, c_min = 0) {
    rel_x <- x / time_len
    {sin({(4 * phases)*rel_x + 1} * {pi/2}) + 1} * 0.5 * (c_max - c_min) + c_min
}
u_curve <- function(time_len, phases = 1, c_max = 1, c_min = 0) {
    sapply(1:time_len, u_curve_1x, time_len, phases, c_max, c_min)
}





# Change a covariance matrix
change_cov_mat <- function(cov_matrix_0, i1, i2, new_cov){
    if (i1 == i2) { stop("Cannot change covariance with itself.") }
    new_cov_matrix <- cov_matrix_0
    new_cov_matrix[i1, i2] <- new_cov
    new_cov_matrix[i2, i1] <- new_cov
    return(new_cov_matrix)
}







# 3 species
# 1:2 decline with X
# 3 increases with X
# when X reaches a threshold, r for 2 affects r for 3, and vice versa


sym_3sp <- function(b2_change_fun, phases = 1, abundances_0 = rep(1000, 3), 
                    b1s = c(-0.1,0.1), phase_len = 100, r_sd = 0.05, K = 2000,
                    spp_names = c('naive', 'host', 'symbiont')) {
    time_len <- as.integer(phases * phase_len)
    X <- u_curve(time_len, phases = phases, c_max = 1, c_min = -1)
    n_spp <- 3
    pop_sizes <- matrix(integer(time_len*n_spp), ncol = n_spp)
    
    pop_sizes[1,] <- abundances_0  # Initial abundances
    
    for (t in 2:time_len) { # t = 2; # t = phase_len/2
        # How sp 3 affects 2, and vice versa
        b2s <- b2_change_fun(t/phase_len)
        # Changing rs_t based on X[t], plus adding random noise
        # Linear deterministic portion based on X[t] and other sp abundances
        lin_det <- b1s[c(1,1,2)] * X[t] + c(0, b2s) * (pop_sizes[t-1,]/ K)
        # Note: `- 0.5` is to make growth rates center around zero
        rs_t <- inv_logit(lin_det) - 0.5 + rnorm(n_spp, sd = r_sd)
        pop_sizes_t <- log_growth(pop_sizes[t-1,], rs_t, K)
        pop_sizes[t,] <- pop_sizes_t
    }
    pop_sizes_df <- as_data_frame(pop_sizes)
    colnames(pop_sizes_df) <- spp_names
    pop_sizes_df <- pop_sizes_df %>%
        mutate(generation = 1:time_len) %>%
        select(generation, everything())
    return(pop_sizes_df)
}





b2_change_fun <- function(rel_t) {
    i <- stepfun(c(0.25, 0.5, 0.75), 1:4)(rel_t)
    # Effect of species 3 on 2
    b3_2 <- c(0, 0.1, 0.5, 0.2)[i]
    # Effect of species 2 on 3
    b2_3 <- c(0, -0.1, 0.1, 0.2)[i]
    return(c(b3_2, b2_3))
}


# Setting simulation parameters
n_phases <- 4
n_sims <- 100
phase_len <- 100
b1s <- c(-0.25, 0.5)
# sim_df <- sym_3sp(b2_change_fun, phases = n_phases, b1s = c(-0.25, 0.5)) %>%
#     gather(species, abundance, -generation, factor_key = TRUE)

# 100 simulations took 2.59 secs
s0 = Sys.time()
set.seed(8)
sim_df <- lapply(1:n_sims, 
                 function(i) {
                     sym_3sp(b2_change_fun, phases = n_phases, b1s = b1s) 
                 }) %>%
    bind_rows %>%
    mutate(sim = rep(1:n_sims, each = phase_len*n_phases)) %>%
    gather(species, abundance, -generation, -sim, factor_key = TRUE)
s1 = Sys.time()
cat(n_sims, 'simulations took', round(s1 - s0, 2) %>% as.numeric, 'secs')



sim_df %>% 
    ggplot(aes(generation, abundance, group = sim)) +
    stat_function(fun = u_curve_1x,
                  args = list(phases = n_phases, time_len = n_phases * 100, 
                              c_max = max(sim_df$abundance),
                              c_min = min(sim_df$abundance)),
                  linetype = 3, size = 0.25) +
    geom_line(alpha = 0.3, color = 'dodgerblue') +
    theme_lan() +
    facet_grid( ~ species)


# Focusing on 1 simulation
sim_df %>% 
    filter(sim == 1) %>% 
    ggplot(aes(generation, abundance, color = species)) +
    stat_function(fun = u_curve_1x,
                  args = list(phases = n_phases, time_len = n_phases * 100, 
                              c_max = max(sim_df$abundance),
                              c_min = min(sim_df$abundance)),
                  linetype = 3) +
    geom_line() +
    theme_lan()



































# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================
# =======================================================================================






































# 
# sim_plot_og <- function(y_mat, phases = 1, phase_len = 1000,
#                      sp_names = c('incr D', 'incr D, host', 'symbiont')){
#     time_len <- as.integer(phase_len * phases)
#     as_data_frame(y_mat) %>% 
#         mutate(generation = seq(time_len)) %>% 
#         gather(species, abundance, -one_of(c('generation'))) %>% 
#         mutate(species = sp_names[gsub('V', '', species) %>% as.integer] %>% factor) %>% 
#         ggplot(aes(generation, abundance, color = species)) + 
#         geom_line(
#             data = data_frame(
#                 generation = seq(time_len),
#                 abundance = u_curve(time_len, phases, 
#                                     c_max = max(y_mat), c_min = min(y_mat)), 
#                 species = unique(factor(sp_names))[1]
#             ), linetype = 3, color = 'black') +
#         geom_line() +
#         theme_lan()
# }
# 
# 
# 
# # Simulates abundance at time t+1
# abun_fun <- function(N_t, p_birth, p_death, cov_matrix = NULL, rando = TRUE) {
#     
#     if (rando) {
#         # (Assuming some dead indiv bred)
#         births <- mapply(function(N_t, p_b){rpois(1, N_t * p_b)}, N_t, p_birth)
#         births <- matrix(births, ncol = 1)
#         deaths <- mapply(function(N_t, p_d){rpois(1, N_t * p_d)}, N_t, p_death)
#         deaths <- matrix(deaths, ncol = 1)
#     } else {
#         # (Assuming some dead indiv bred)
#         births <- N_t * p_birth
#         births <- matrix(births, ncol = 1)
#         deaths <- N_t * p_death
#         deaths <- matrix(deaths, ncol = 1)
#     }
#     
#     nvals <- length(N_t)
#     
#     # Including covariance of death rates
#     if (!is.null(cov_matrix)) {
#         cov_dc <- chol(cov_matrix)  # Cholesky decomposition of cov matrix
#         new_deaths_1 <- t(t(cov_dc) %*% deaths)  # <-- Actually including cov
#         new_deaths_2 <- t(t(cov_dc) %*% deaths[nvals:1,])[,nvals:1]
#         new_deaths_2 <- matrix(new_deaths_2, nrow = 1)
#         new_deaths <- apply(array(c(p,q), dim = c(dim(new_deaths_1),2)), c(1,3), min)
#         # new_deaths <- (new_deaths_1 + new_deaths_2)/2
#     } else {
#         new_deaths <- deaths
#     }
#     y_t <- sapply(N_t - deaths + births, function(x){max(c(1,x))})
#     y_t <- matrix(y_t, nrow = 1)
#     return(y_t)
# }
# 
# 
# 
# 
# 
# 
# # 3 species
# # 1:2 decline with X
# # 3 doesn't depend on it
# # when X reaches a threshold, p_death for 2 covaries with 3
# 
# 
# sym_sim <- function(phases, abundances_0 = rep(1000, 3), b1s = c(1,1), 
#                     phase_len = 1000L, p_birth = 0.25, covs = c(0, 0.25, 0.75)) {
#     
#     time_len <- as.integer(phases * phase_len)
#     X <- u_curve(time_len, phases = phases, c_max = 1, c_min = 0)
#     n_spp <- length(abundances_0)
#     y <- matrix(integer(time_len*n_spp), ncol = n_spp)
#     
#     y[1,] <- abundances_0  # Initial abundance
#     p_death_even <- logit(p_birth)
#     
#     # Covariance between death rates to start with
#     cov_mat_0 <- matrix(c(1.0,    0.0,    0.0,
#                           0.0,    1.0,    0.0,
#                           0.0,    0.0,    1.0),
#                         ncol = n_spp)
#     
#     p_deaths <- rep(p_death_even, 3)
#     cov_matrix <- cov_mat_0
#     for (t in 2:time_len) {
#         # Covariance inversely related to X[t] bc at high X, death low so no selection
#         # for symbiosis
#         cov_t <- stepfun(time_len/phases * c(0.25, 0.5), covs)(t)
#         cov_matrix <- change_cov_mat(cov_matrix, 2, 3, cov_t)
#         # High X reduces death rates
#         p_deaths[1:2] <- inv_logit(p_death_even + (b1s[1] * X[t] / 1e3))
#         p_deaths[3] <- inv_logit(p_death_even - (b1s[2] * X[t] / 1e3))
#         y_t <- abun_fun(y[t-1,], p_birth, p_deaths, cov_matrix = cov_matrix)
#         y[t,] <- y_t
#     }
#     
#     return(as_data_frame(y))
# }
# 
# 
# 
# 
# 
# 
# n_phases <- 0.5
# n_sims <- 100
# phase_len <- 2e2L
# 
# 
# df <- lapply(1:n_sims, 
#              function(i) sym_sim(n_phases, b1s = c(1,2), phase_len = phase_len, 
#                                  covs = c(0, 0.9, 0.9999))) %>% 
#     bind_rows %>% 
#     mutate(sim = rep(1:n_sims, each = phase_len*n_phases),
#            generation = rep(1:(phase_len*n_phases), n_sims)) %>% 
#     gather(species, abundance, -one_of(c('generation', 'sim'))) %>% 
#     mutate(species = sp_names[gsub('V', '', species) %>% as.integer] %>% factor)
# 
# 
# df %>% 
#     ggplot(aes(generation, abundance, group = sim)) + 
#     # stat_function(fun = u_curve_p, 
#     #               args = list(n_phases*phase_len, n_phases, 
#     #                           c_max = max(df$abundance), c_min = min(df$abundance)),
#     #               linetype = 3) +
#     geom_line(alpha = 0.3, color = 'dodgerblue') +
#     theme_lan() +
#     facet_grid( ~ species)
# 
# 
#     
# 
# sim_plot <- function(y_df, phases = 1, phase_len = 1000,
#                         sp_names = c('incr D', 'incr D, host', 'symbiont')){
#     time_len <- as.integer(phase_len * phases)
#     as_data_frame(y_mat) %>% 
#         mutate(generation = seq(time_len)) %>% 
#         gather(species, abundance, -one_of(c('generation'))) %>% 
#         mutate(species = sp_names[gsub('V', '', species) %>% as.integer] %>% factor) %>% 
#         ggplot(aes(generation, abundance, color = species)) + 
#         geom_line(
#             data = data_frame(
#                 generation = seq(time_len),
#                 abundance = u_curve(time_len, phases, 
#                                     c_max = max(y_mat), c_min = min(y_mat)), 
#                 species = unique(factor(sp_names))[1]
#             ), linetype = 3, color = 'black') +
#         geom_line() +
#         theme_lan()
# }



# y %>% as_data_frame

# inv_logit(b1 * X[1])














# 
# 
# 
# # ~~~~~~~~~~~~~~~~~~~
# # The following (inside "~~~") was from 
# # https://www.r-bloggers.com/simulating-data-following-a-given-covariance-structure/
# 
# # Between genes in species 1
# covM1_0 <- matrix(c(1.0,    0.0,    0.0,    0.0,
#                     0.0,    1.0,    0.0,    0.0,
#                     0.0,    0.0,    1.0,    0.0,
#                     0.0,    0.0,    0.0,    1.0),
#                   nrow = 4, ncol = 4) %>% 
#     change_cov_mat(., 1, 2, 0.5) %>% 
#     change_cov_mat(., 1, 3, 0.25) %>% 
#     change_cov_mat(., 2, 3, 0.1) %>% 
#     change_cov_mat(., 2, 4, 0.75)
# # Between genes in species 2
# covM2_0 <- matrix(c(1.0,    0.0,    0.0,    0.0,
#                     0.0,    1.0,    0.0,    0.0,
#                     0.0,    0.0,    1.0,    0.0,
#                     0.0,    0.0,    0.0,    1.0),
#                   nrow = 4, ncol = 4) %>% 
#     change_cov_mat(., 1, 2, 0.5) %>% 
#     change_cov_mat(., 1, 3, 0.25) %>% 
#     change_cov_mat(., 2, 3, 0.1) %>% 
#     change_cov_mat(., 2, 4, 0.75)
# 
# 
# # Between genes in all species
# covM_0 <- matrix(c(
#     #           sp1              ||         sp2
#     1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
#     0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
#     0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,
#     0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,
#     # -----                      ||
#     0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,
#     0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,
#     0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,
#     0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0
#     ), nrow = 8, ncol = 8)
# 
# # # colnames(covM_0) <- c("sp1_g1", "sp1_g2", "sp1_g3", "sp1_g4", 
# # #                       "sp2_g1", "sp2_g2", "sp2_g3", "sp2_g4")
# # 
# 
# log_growth <- function(x, K = 1, r = 0.1, x_0 = 50) {
#     K / {1 + exp(-r * (x - x_0))}
# }
# 
# int_log_growth <- function(Nt, r, K = 100) {
#     Nt + {r*Nt * (1 - Nt/K)}
# }
# 
# 
# 
# time_len = 100
# 
# X <- seq(time_len)
# 
# # Net fitness contributions
# b1 <- c(0.05, 0.1, 0.25, 0.6)
# a1 = 0.5
# 
# 
# nvars <- dim(covM_0)[1]
# 
# 
# 
# covM <- covM_0
# covM1 <- covM1_0
# covM2 <- covM2_0
# 
# y <- matrix(numeric(nvars * time_len), ncol = nvars)
# y[1,] <- 1
# 
# for (t in 2:time_len){
#     if (t == as.integer(time_len/2)) {
#         covM <- change_cov_mat(covM, 1, 5, 0.75)
#     }
#     # Cholesky decomposition
#     covL1 <- chol(covM1)
#     # Without covariance...
#     # y1_t <- a1 * y[t-1,1:4] + b1 + rnorm(nvars/2)
#     y1_t <- a1[1] * y[t-1,1:4] * b1 + rnorm(nvars/2)
#     y1_t <- matrix(y1_t, nrow = nvars/2)
#     # Including covariance
#     y1_t <- t(t(covL1) %*% y1_t)
#     
#     # Cholesky decomposition
#     covL2 <- chol(covM2)
#     # Without covariance...
#     # y2_t <- a1 * y[t-1,5:8] + b1 + rnorm(nvars/2)
#     y2_t <- a1[2] * y[t-1,5:8] * b1 + rnorm(nvars/2)
#     y2_t <- matrix(y2_t, nrow = nvars/2)
#     # Including covariance
#     y2_t <- t(t(covL2) %*% y2_t)
#     
#     # Cholesky decomposition
#     covL <- chol(covM)
#     # Including covariance
#     y[t,] <- t(t(covL) %*% c(y1_t, y2_t))
# }; rm(t)
# 
# # int_log_growth(0.5, r = 0.7)
# 
# 
# 
# # for (t in 2:time_len){
# #     covM <- change_cov_mat(covM, 1, 4, inv_logit(5 * t/time_len))
# #     covM <- change_cov_mat(covM, 5, 8, inv_logit(2 * t/time_len))
# #     covM <- change_cov_mat(covM, 2, 6, inv_logit(5 * t/time_len))
# #     # Cholesky decomposition
# #     covL <- chol(covM)
# #     # Without covariance...
# #     y_t <- a1 * y[t-1,] + b1 * X[t] + rnorm(nvars)  # MASS::mvrnorm(1, rep(0,8), covM)
# #     y_t <- matrix(y_t, nrow = nvars)
# #     
# #     # Including covariance
# #     y[t,] <- t(t(covL) %*% y_t)
# # }
# 
# colnames(y) <- c("sp1_g1", "sp1_g2", "sp1_g3", "sp1_g4", 
#                  "sp2_g1", "sp2_g2", "sp2_g3", "sp2_g4")
# df <- as_data_frame(y) %>% 
#     mutate(generation = seq(time_len)) %>% 
#     gather(sp_gene, abundance, -one_of(c('generation'))) %>% 
#     mutate(
#         species = factor(ifelse(grepl('sp1', sp_gene), 'species 1', 'species 2')),
#         gene = factor(substr(sp_gene, 6, 6))
#         )
# 
# df %>% 
#     group_by(species, generation) %>% 
#     summarize(abundance = sum(abundance)) %>% 
#     ggplot(aes(generation, abundance)) +
#     geom_line() +
#     facet_grid(~species) +
#     theme_lan()
# 
# 
# 
# 
# 
# # ~~~~~~~~~~~~~~~~~~~
# 
# for (t in 1:nsims){
#     X[t] <- a1 * X[t-1] + b1 * U[t] + rnorm(1)
# }