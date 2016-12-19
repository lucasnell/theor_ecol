
library(tidyverse)
library(lme4)
library(lmerTest)
library(GGally)
library(vegan)

# Custom ggplot2 theme
theme_lan <- function(base_size = 10, base_family = 'Helvetica') {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'bold', size = 11),
            strip.background = element_rect(color = NA, fill = NA),
            panel.grid = element_blank()
        )
}

z_trans <- function(x) {(x - mean(x)) / sd(x)}
min_log <- function(x) {
    nz <- min(x[x>0])
    y <- log(x + nz)
    return(y)
}
anti_min_log <- function(y, x_orig){
    nz <- min(x_orig[x_orig>0])
    x <- exp(y) - nz
    return(x)
}
asin_sqrt <- function(x){
    y <- asin(sqrt(x))
    return(y)
}
anti_asin_sqrt <- function(y){
    x <- sin(y)^2
    return(x)
}



# This problem set is based on the Pine Barrens plant comunity data
# in Li & Waller (2015); a pdf of the paper is in the "Pine barrens data" 
# folder. The data set consists of 30 sites sampled in 2012. At each site,
# 50 quadrats were sampled, with presence/absence of species scored. 
# Therefore, species abundances are given on a scale of counts from 0 to 50.

## Environmental variables using `rda`

##  Import data "veg_aggregated_1958.csv" and "veg_aggregated_2012.csv"
veg_df <- read_csv(paste0('./box/ZooEnt_540_2016/Data set analyses/',
                          'Pine barrens data 28Nov16/veg_aggregated_2012.csv'),
                   col_types = 'cci')


# d = veg_df -->
# w = veg_df_wide -->


# Reducing the number of species.
sp_to_keep <- {veg_df %>% 
    group_by(sp) %>% 
    summarize(occurs = sum(count > 0), mean_count = mean(count)) %>% 
    filter(occurs >= 5, mean_count >= 0.5)}$sp

veg_df <- veg_df %>% 
    filter(sp %in% sp_to_keep) %>% 
    mutate(den = asin_sqrt(count/50))


##  Import data "envi.csv"
env_df <- read_csv(paste0('./box/ZooEnt_540_2016/Data set analyses/',
                     'Pine barrens data 28Nov16/envi.csv'),
                   col_types = paste(c('ci', rep('d', 30), 'c'), collapse = '')) %>% 
    filter(date == 2012)


# Add a subset of environmental variables to d; you can add more variables if you want
veg_df$tmin <- env_df$tmin[match(veg_df$site, env_df$site)]
veg_df$H <- env_df$H.prop[match(veg_df$site, env_df$site)]
veg_df$Mn <- env_df$Mn.mg_kg[match(veg_df$site, env_df$site)]
veg_df$shade <- env_df$shade[match(veg_df$site, env_df$site)]
veg_df$precip <- env_df$precip[match(veg_df$site, env_df$site)]

# Standardize environmental variables to have mean = 0 and variance = 1 so that the
# coefficients can be compared. This also makes analyses more stable (less likely to 
# have convergence problems).

veg_df <- veg_df %>% 
    mutate_at(vars(tmin, H, Mn, shade, precip), funs(z_trans)) %>% 
    mutate_at(vars(site, sp), funs(factor))

veg_df_wide <- veg_df %>% select(-count) %>% spread(sp, den)



#########################################################################
# Question 1: Analyze data using RDA. Which environmental variables seem
# to be important for determining the composition of the communities? 
# To answer this question, you can use the code from the pine beetle
# analyses, remembering that the environmental variables are much
# easier to handle than "year" in the pine beetle data set since
# they are structurally independent.
#########################################################################



nsite <- nlevels(veg_df$site)
nspp <- nlevels(veg_df$sp)

# create community matrix
comm_matrix <- matrix(veg_df$den, c(nsite, nspp), byrow = TRUE)
colnames(comm_matrix) <- levels(veg_df$sp)
rownames(comm_matrix) <- 1:nsite


envir_matrix <- cbind(veg_df_wide$tmin, veg_df_wide$H, veg_df_wide$Mn, 
                      veg_df_wide$shade, veg_df_wide$precip)
colnames(envir_matrix) <- c('tmin', 'H', 'Mn', 'shade', 'precip')
rownames(envir_matrix) <- veg_df_wide$site



# perform PCA
pca_fit <- rda(comm_matrix, scale = TRUE)
# plot(pca_fit, type="text")

# fit rda
rda_fit <- rda(comm_matrix, envir_matrix, scale = TRUE)
plot(rda_fit, type="text")


# ggpairs(veg_df_wide, columns = 2:6)


# get partial p-values that give differences of years from year2004
set.seed(999)
rda_partial_pvals <- sapply(1:ncol(envir_matrix), function(j){
    # partial rda
    Z <- envir_matrix[,-j]
    rda_fit_partial <- rda(comm_matrix, envir_matrix[,j], Z, scale = TRUE)
    anova_rda <- permutest(rda_fit_partial, permutations=10000)
    pval <- mean(anova_rda$F.0 < anova_rda$F.perm)
    pval
})
names(rda_partial_pvals) <- colnames(envir_matrix)
# tmin      H     Mn  shade precip 
# 0.0238 0.8114 0.3381 0.3540 0.0096






## Environmental variables using `lmer`



############################################################################
# Question 2: Analyze the community data using lmer (multilevel model) 
# to test whether the environmental variables affect community composition.
# One thing to consider is whether you test each environmental variable
# separately (univariate analyses) or together (multivariate analyses). 
############################################################################


env_vars <- colnames(envir_matrix)

all_mods <- lapply(1:length(env_vars), function(i) {
    var_mat <- t(combn(env_vars, i))
    apply(var_mat, 1, function(x) {
        paste(c('(1 | sp)', as.character(x), paste0('(0 + ', as.character(x), ' | sp)')), 
              collapse = ' + ')
    })
    }) %>% 
    c(., recursive = TRUE)


aics <- sapply(all_mods, 
               function(m) {
                   AIC(lmer(as.formula(paste('den ~', m)), data = veg_df, REML = FALSE))
                }, USE.NAMES = FALSE)

# sort(aics) - min(aics)

# I was going to use the formula for the model with lowest AIC, but just ended up using
# the full model
# final_form <- as.formula(paste('den ~', all_mods[aics == min(aics)]))
final_form <- as.formula(paste('den ~', all_mods[length(all_mods)]))

lmer_mod <- lmer(final_form, data = veg_df)
# summary(lmer_mod)

# hist((veg_df$den))


p_strings <- c('< 0.001', '= 1.0', '= 0.022', '= 0.031', '< 0.001')
cat(paste(env_vars, p_strings, collapse = '\n'))




## Using `lmer` to perform ordination


############################################################################
# Question 3: Use lmer to perform the "ordination" analyses as in
# Jackson et al. (2012). Is the biplot for the MLM similar to the
# RDA biplot? For this, as a template you can use the code from
# Jackson et al. (2012) that is provided in the "Pine barrens" folder.
############################################################################

library(stringr)
no_ran_form <- str_sub(paste(final_form)[3], start = 1L, 
                   end = str_locate(paste(final_form)[3], '0 +')[[1,'start']] - 5) %>% 
    paste('den ~', .)
no_ran_form <- as.formula(no_ran_form)

lmer_mod_nr <- lmer(no_ran_form, data = veg_df)


MLM_fitted <- array(fitted(lmer_mod) - fitted(lmer_mod_nr), c(nsite, nspp))
rownames(MLM_fitted) = c(1:nrow(MLM_fitted))
colnames(MLM_fitted) = levels(veg_df$sp)

# head(MLM_fitted)

# standardize over spp
MLM_fitted_standard <- lapply(1:nspp, function(j) z_trans(MLM_fitted[, j])) %>% 
    do.call(cbind, .)


ss <- cor(MLM_fitted_standard)
U <- svd(ss)
mlm_fit <- MLM_fitted_standard %*% U$v
mlm_fit <- mlm_fit[, 1:2]



# environmental variables (only those with random effects)

env_mat <- veg_df_wide %>% select(one_of(env_vars)) %>% as.matrix

mlm_envir <- NULL
for (j in 1:length(env_vars)){
    mlm_envir <-
    cbind(mlm_envir, env_mat[, j] * mlm_fit[, 1], env_mat[, j] * mlm_fit[, 2])
}
envir_points <- matrix(colMeans(mlm_envir), nrow = ncol(env_mat), ncol = 2, byrow = T)



par(mfrow = c(1, 2))

plot(rda_fit, type="text", main = "RDA")

# plot mlm
plot(-mlm_fit,
     xlab = "PC1",
     ylab = "PC2",
     type = "n",
     main = "lmer")
text(-mlm_fit, label = 1:nsite, cex = 0.5)

arrow_coordMLM <- cbind(array(0, dim(envir_points)), -envir_points)

arrows(
    arrow_coordMLM[, 1],
    arrow_coordMLM[, 2],
    arrow_coordMLM[, 3],
    arrow_coordMLM[, 4],
    col = "black",
    length = 0.1
)

text(
    1.3 * -envir_points,
    label = env_vars,
    cex = 0.7
)

mlm_scores <- mlm_fit
rda_scores <- scores(rda_fit, display = "lc")

proc <- procrustes(mlm_scores, rda_scores,
                   scale = TRUE, symmetric = TRUE)



## Procrustes sum of squares: 
signif(proc$ss, 3)




# Moving this file to Box folder...
system(
    paste("cd", getwd(),
          "&& cp Nell_pine_barrens.R",
          "~/'Box Sync/ZooEnt_540_2016/Data set analyses/Pine barrens data 28Nov16'")
)
