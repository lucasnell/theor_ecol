####  Meta Data ####

# Data collected between 2004-2008 in 24 sites in 4 regions in Wisconsin
# CS =  Central Sands
# KM =  Kettle Morraine
# SG =  Spring Green
# WS =  West Salem

# Columns 10:27 are various insect species. I_pini and I_grandicollis are of particular interest. 


# Traps sampled every two weeks from May-September

# 3 Treatments: Root Sever, Pocket Control, Asymptomatic
# 3 Trap Types: Funnel, Milk Jug, and Pitfall Traps

## Clear
# rm(list=ls())

####  Packages  ####
library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)
# library(nlme)

z_trans <- function(x){
    (x - mean(x)) / sd(x)
}

log_min <- function(x){
    nz <- min(x[x > 0])
    log(x + nz)
}

min_trans <- function(x){
    x_min <- min(x)
    x - x_min
}


####  Read in Data  ####


reg_sites <- bind_rows(
    expand.grid('CS', c('bal', 'bkr', 'blkr', 'sch')),
    expand.grid('KM', c('ktm')),
    expand.grid('SG', c('bp', 'gat', 'scf')),
    expand.grid('WS', c('wsm')))


raw_dat <- read_csv(paste0('~/GitHub/Wisconsin/theor_ecol/box/ZooEnt_540_2016/',
                          'Data set analyses/Red Pine Decline data 14Nov16/', 
                          'rpd_data.csv')) %>% 
    mutate(date = as.Date(date, format = '%m/%d/%y'),
           region = sapply(site, function(s){reg_sites[reg_sites[,2]== s, 1]})) %>% 
    gather(species, count, i_pini:p_pallens)



# rpd_df <- raw_dat %>% 
#     select(-lure, -trap_type) %>% 
#     group_by(date, site, subsite, treatment, region, species) %>% 
#     summarize(
#         year = as.integer(median(year)),
#         count = sum(count)
#     ) %>% 
#     ungroup %>% 
#     filter(!is.na(count))

rpd_dfm <- raw_dat %>% 
    group_by(date, site, subsite, treatment, region, species, trap_type) %>% 
    summarize(
        year = paste(median(year)),
        count = mean(count)
    ) %>% 
    ungroup %>% 
    filter(!is.na(count)) %>% 
    mutate_if(is.character, factor)



# AC = asymptomatic control (healthy pine stands)
# PC = "pocket control" (??)
# RS = root sever treatment


# spc=a[,7:24]
# ord=metaMDS(spc, autotransform = TRUE, trymax=1000, distance = "manhattan")
# plot(ord)
# text(ord, display="spec", col="red")


# Takes ~30 seconds
# lme_mod <- glmer(count ~ (1 | species) + treatment + (treatment | species) + (1 | site), 
#                  family = poisson, data = rpd_df)
# lme_mod_r0 <- glmer(count ~ (1 | species) + treatment, family = poisson, data = rpd_df)
# lme_mod_f0 <- glmer(count ~ (1 | species), family = poisson, data = rpd_df)

# # Below takes about 5 minutes
# lme_mod <- lmer(log_min(count) ~ (1 | species) +  # (1 | subsite) + 
#                     treatment + (treatment | species) + 
#                     # year + (year | species) +
#                     yr_trans + (yr_trans | species) +
#                     trap_type + (trap_type | species),  # + treatment:trap_type,
#                 data = rpd_dfm)
# 
# 
# s0 = Sys.time()
# lme_mod2 <- lmer(log_min(count) ~ (1 | species) +  # (1 | subsite) + 
#                      treatment + (treatment | species) + 
#                      year + (year | species) +
#                      # yr_trans + (yr_trans | species) +
#                      trap_type + (trap_type | species),
#                  control = lmerControl(optimizer = "nloptwrap"),
#                  data = rpd_dfm)
# s1 = Sys.time()
# s1 - s0
# 
# rpd_dfm$year %>% unique
# 
# str(lme_mod2)
# 
# 
# 
# # mlm <- lmer(log(qq$ABUNDANCE+min) ~ (1|SPP)+year+(year|SPP)+trap+(trap|SPP)+
# #                 treatment+(treatment|SPP),
# #             control=lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)), 
# #             data=qq)
# 
# 
# # AIC(lme_mod)
# 
# 
# # P-value for including random effect of treatment
# pchisq(2*(as.numeric(logLik(lme_mod)) - as.numeric(logLik(lme_mod_r0))), df=1,
#        lower.tail = FALSE)
# # P-value for including fixed effect of treatment
# pchisq(2*(as.numeric(logLik(lme_mod_r0)) - as.numeric(logLik(lme_mod_f0))), df=1,
#        lower.tail = FALSE)
# 
# 
# summary(lme_mod)
# 
# coef(lme_mod)$species
# 
# ranef(lme_mod)$species
# 
# fixef(lme_mod)
# 
# 
# 
# plot(lme_mod)
# 
# data_frame(residuals = resid(lme_mod), x = factor(rpd_dfm$treatment)) %>% 
#     ggplot(aes(x, residuals)) +
#     geom_jitter(alpha = 0.1, shape = 1) +
#     stat_summary(fun.data = "mean_cl_boot", size = 0.375, color = 'red',
#                  geom = 'errorbar') +
#     theme_bw()
# 
# 
# expand.grid(treatment = unique(rpd_df$treatment), 
#                   species = unique(rpd_df$species), stringsAsFactors = FALSE) %>% 
#     as.tbl %>% 
#     mutate(count = predict(lme_mod, newdata = nd)) %>% 
#     ggplot(aes(as.numeric(factor(treatment)), count, color = species)) + 
#     geom_line()



'log_min(count) ~'
mod_vars <- lapply(1:3, function(i) {
    mat <- t(combn(c('treatment', 'year', 'trap_type'), i))
    out <- character(nrow(mat))
    for (j in 1:nrow(mat)){
        out[j] <- paste(mat[j,], collapse = ' + ')
    }
    return(out)
    }) %>% 
    c(., recursive = TRUE)
mod_vars <- c(mod_vars, 'trap_type * treatment', 'year + trap_type * treatment')

one_sp_mod <- function(sp) {
    df <- rpd_dfm %>% filter(species == sp)
    mods <- lapply(mod_vars,
                   function(m) {
                       form <- as.formula(paste('log_min(count) ~ (1 | subsite) + ', m))
                       lmer(form, data = df, REML = FALSE)
                   })
    mod_aics <- data_frame(mod = mod_vars, aic = sapply(mods, AIC), species = sp)
    return(mod_aics)
}

aics <- lapply(levels(rpd_dfm$species), one_sp_mod) %>% 
    bind_rows


aics %>% group_by(mod) %>% summarize(aic = median(aic)) %>% 
    arrange(aic) %>% mutate(d_aic = aic - min(aic))


# # This takes a while but is about how it's done
# lmer(log_min(count) ~ (1 | subsite) + year + treatment * trap_type +
#          (0 + year | species) + (0 + trap_type | species) + (0 + treatment | species) + 
#          (0 + trap_type + treatment | species), 
#      # control = lmerControl(optimizer = "nloptwrap"),
#      data = rpd_dfm)


# 
# # Moving this file to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp Nell_red_pine.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/'")
# )