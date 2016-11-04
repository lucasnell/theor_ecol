library(magrittr)
# library(ggthemes)
library(tidyverse)
library(lme4)
library(car)
library(nlme)
library(mgcv)
library(broom)

if (grepl('theor_ecol$', getwd())){
    setwd('./Pine_dat')
} else if (!grepl('Pine_dat$', getwd())) {
    stop('Unknown directory. It should be in one ending in "theor_ecol" or "Pine_dat".')
}
# Custom ggplot2 theme
theme_lan <- function(base_size = 10, base_family = 'Helvetica') {
    ggthemes::theme_fivethirtyeight(base_size = base_size, 
                                    base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'bold'),
            axis.title = element_text(size = 11)
        )
}


# obligatory floating egg stage --> lots of mixing, perhaps big floods are bad


catch_df <- read_csv('catch_data.csv', col_types = 'iccdidi') %>% 
    mutate(area = ifelse(area == 0, NA, area)) %>% 
    filter(isolated_pools == 0) %>% 
    select(-isolated_pools)






z_trans <- function(vec, pre_trans = NULL){
    if (!is.null(pre_trans)){
        if (is.character(pre_trans)) {
            new_vec <- eval(as.name(pre_trans))(vec)
        } else if (is.function(pre_trans)) {
            new_vec <- pre_trans(vec)
        } else {
            stop('pre_trans must be a character or function.')
        }
    } else {
        new_vec <- vec
    }
    (new_vec - mean(new_vec)) / sd(new_vec)
}


log_min <- function(vec, type = log) {
    new_vec <- vec + min(vec[vec > 0])
    if (is.character(type)) {
        logged_vec <- eval(as.name(type))(new_vec)
    } else if (is.function(type)) {
        logged_vec <- type(new_vec)
    } else {
        stop('type must be a character or function.')
    }
    return(logged_vec)
}


flow_df <- read_csv('all_flow_data.csv', 
                    col_types = 'icddiiiiii') %>% 
    mutate(z_mean = z_trans(mean, pre_trans = log)) %>% 
    gather(key = threshold, value = duration, starts_with('over_')) %>% 
    mutate(threshold = as.integer(gsub('over_', '', threshold))) %>% 
    group_by(threshold) %>% 
    mutate(z_duration = z_trans(duration, pre_trans = log_min)) %>% 
    ungroup



# Only meant for one threshold
get_flows <- function(year, reach, get_z = FALSE, thresh = 3000){
    column_name <- ifelse(get_z, 'z_duration', 'duration')
    one_flow <- function(y, r, t, c){
        flow <- as.numeric(flow_df[flow_df$year == y & flow_df$reach == r & 
                                       flow_df$threshold == t, ][[c]])
        ifelse(length(flow) == 0, NA, flow)
    }
    flows <- mapply(one_flow, year, reach, t = thresh, c = column_name)
    return(flows)
}

get_mean <- function(year, reach, get_z = FALSE){
    column_name <- ifelse(get_z, 'z_mean', 'mean')
    one_flow <- function(y, r, c){
        flow <- as.numeric(flow_df[flow_df$year == y & flow_df$reach == r & 
                                       flow_df$threshold == 1000, ][[c]])
        ifelse(length(flow) == 0, NA, flow)
    }
    flows <- mapply(one_flow, year, reach, c = column_name)
    return(flows)
}


# Aggregate by reach and year
yr_reach <- catch_df %>% 
    group_by(year, reach) %>% 
    summarize(dry = ifelse(any(is.na(area)), 1L, 0L),
              count = sum(count, na.rm = TRUE), 
              area = sum(area, na.rm = TRUE)) %>% 
    ungroup %>%
    mutate(flow = get_flows(year, reach, get_z = FALSE),
           flow_z = get_flows(year, reach, get_z = TRUE),
           mean_z = get_flows(year, reach, get_z = TRUE),
           reach = as.factor(reach),
           l.area = log(area),
           id = paste(reach, year, sep = '_'),
           cpue = count / (area/100)) %>%
    arrange(reach, year) %>%
    filter(!is.na(flow))



yr_reach %>% 
    group_by(reach) %>% 
    mutate(count.stan = count/area / max(count/area), flow.stan = flow / max(flow)) %>% 
    ungroup %>% 
    ggplot(aes(year, count.stan)) +
    geom_line() +
    geom_line(aes(y = flow.stan), color = 'red', linetype = 2) +
    facet_grid(reach ~ ., scales = 'free_y') +
    theme_lan() +
    ylab('Standardized CPUE (black) and flow (red)')



catch_df %>% 
    group_by(year) %>% 
    summarize(count = sum(count, na.rm = TRUE), 
              area = sum(area, na.rm = TRUE),
              flow = get_flows(year, reach = 'Angostura', get_z = FALSE)[1]) %>% 
    ungroup %>%
    mutate(flow.stan = flow / max(flow), count.stan = (count/area) / max(count/area)) %>% 
    ggplot(aes(year, count.stan)) +
    geom_line() +
    geom_line(aes(y = flow.stan), color = 'red', linetype = 2) +
    theme_lan() +
    ylab('Standardized CPUE (black) and flow (red)')










rm_df <- data_frame(starts = catch_df$rm_start %>% unique,
                    appears = sapply(
                        starts, 
                        function(s){nrow(catch_df[catch_df$rm_start == s,])}
                        )) %>% 
    arrange(desc(appears))
rm_df %>% as.data.frame

starts <- rm_df$starts[rm_df$appears >= 6]

yr_reach_filt <- catch_df %>% 
    filter(rm_start %in% starts) %>% 
    group_by(year, reach) %>% 
    summarize(dry = ifelse(any(is.na(area)), 1L, 0L),
              count = sum(count, na.rm = TRUE), 
              area = sum(area, na.rm = TRUE)) %>% 
    ungroup %>%
    mutate(flow = get_flows(year, reach, get_z = TRUE),
           reach = as.factor(reach),
           l.area = log(area),
           id = paste(reach, year, sep = '_')) %>%
    arrange(reach, year) %>%
    filter(!is.na(flow))








yr_reach %>% 
    ggplot(aes(flow, log_min(count/area))) + 
    geom_point(aes(color = factor(reach))) +
    stat_smooth(method = 'glm', formula = y ~ x + I(x^2), se = FALSE) +
    theme_lan()

catch_df %>% 
    group_by(year, reach) %>% 
    summarize_at(vars(count, area), funs(sum), na.rm = TRUE) %>% 
    ungroup %>% 
    mutate(flow = get_flows(year, reach, thresh = 1000, get_z = TRUE)) %>% 
    ggplot(aes((flow), log_min(count/area), color = factor(reach))) + geom_point() +
    theme_lan()


# library(lmerTest)
# lmm <- lmer(count ~ flow + (flow | reach) + offset(area), data = yr_reach)
# summary(lmm)
# Anova(lmm)



# Accounting for temporal autocorrelation
lme_ar1 <- lme(log_min(count) ~ offset(l.area) + flow + dry, 
               random = ~ 1 | reach,  
               correlation = corARMA(form = ~ 1 | reach, p = 1), 
               data = yr_reach, method = 'ML')
AIC(lme_ar1)
plot(lme_ar1)

plot(ACF(lme_ar1, resType = 'normalized'), alpha = 0.5)
# Despite above plot, AIC indicates you should use ar1


# -----------------
# Worse fitting models: no AR and AR2
# -----------------
# lme_mod <- lme(log_min(count) ~ offset(l.area) + flow + I(flow^2),
#                random = ~ 1|reach,  data = yr_reach, method = 'ML')
# AIC(lme_mod)
# # [1] 283.2809
# lme_ar2 <- lme(log_min(count) ~ offset(l.area) + flow + I(flow^2),
#                random = ~ 1|reach,  correlation = corARMA(form = ~ 1 | reach, p = 2),
#                data = yr_reach, method = 'ML')
# AIC(lme_ar2)
# # [1] 258.221


gamm_mod <- gamm(log_min(count) ~ offset(l.area) + s(flow), 
               random = list(reach = ~ 1),  
               correlation = corARMA(form = ~ 1 | reach, p = 1), 
               data = yr_reach, method = 'ML')
summary(gamm_mod$gam)
AIC(gamm_mod$lme)
AIC(lme_ar1)

plot(gamm_mod$lme)
plot(ACF(gamm_mod$lme, resType = 'normalized'), alpha = 0.5)
plot(ACF(lme_ar1, resType = 'normalized'), alpha = 0.5)


# Poisson incorporating a random term that's unique for each observation (`id` column)
gamm_mod_g <- gamm(count ~ offset(l.area) + s(flow), 
                 random = list(reach=~1, id=~1),  
                 correlation = corARMA(form = ~ 1 | reach, p = 1), 
                 data = yr_reach, method = 'ML', family = poisson,
                 niterPQL = 100)

summary(gamm_mod_g$gam)

plot(gamm_mod_g$gam)

nd <- expand.grid(flow_z = seq(-1, 1.5, length.out = 100),
                 l.area = mean(yr_reach$l.area),
                 area = mean(yr_reach$area),
                 reach = yr_reach$reach %>% unique,
                 dry = c(0,1)) %>% 
    as.tbl

# plot(count ~ flow, data = yr_reach)
# lines(predict(gamm_mod_g$gam, newdata = nd, type = 'response') ~ nd$flow)
# lines({exp(predict(gamm_mod$gam, newdata = nd, type = 'response')) - 1} ~ nd$flow,
#       col = 'red')
# 
# 
# 
# plot(gamm_mod$lme)
# plot(gamm_mod_g$lme)



catch_df %>% 
    mutate(flow = get_flows(year, reach, get_z = TRUE),
           cpue = count/area) %>% 
    ggplot(aes(flow, cpue)) + geom_point()



gls_mod <- gls(count ~ offset(area) + flow_z + dry, 
               correlation = corARMA(form = ~ 1 | reach, p = 1), 
               data = yr_reach, method = 'ML', weights = varPower())
AIC(gls_mod)
summary(gls_mod)

data_frame(resid = resid(gls_mod)[order(yr_reach$year)], 
           year = yr_reach$year[order(yr_reach$year)],
           reach = yr_reach$reach[order(yr_reach$year)]) %>% 
    ggplot(aes(year, resid, color = reach)) + 
    geom_point() + geom_line()





AIC(gls(count ~ offset(area) + mean_z + dry, 
        correlation = corARMA(form = ~ 1 | reach, p = 1), 
        data = yr_reach, method = 'ML', weights = varPower()))


yr_reach %>% ggplot(aes(mean_z, flow_z)) + geom_point()

gls_mod2 <- gls(cpue ~ flow_z + dry, 
               correlation = corARMA(form = ~ 1 | reach, p = 1),
               data = yr_reach, method = 'REML', weights = varPower())
AIC(gls_mod2)
summary(gls_mod2)


plot(gls_mod)

plot(gls_mod2)


nd$count <- predict(gls_mod, newdata = nd, type = 'response')
nd$cpue <- predict(gls_mod2, newdata = nd, type = 'response')



ggplot(nd, aes(flow_z, count/area, color = factor(dry))) + 
    geom_point(aes(shape = reach), data = yr_reach) + 
    geom_line() +
    theme_lan()

ggplot(nd, aes(flow_z, cpue, linetype = factor(dry))) + 
    geom_point(aes(color = factor(year)), data = yr_reach) + 
    geom_line() +
    theme_lan()


yr_reach %>% filter(dry == 1)



yr_reach %>% 
    filter(flow_z < -0.5, cpue > 10) %>% 
    ggplot(aes(flow_z, cpue, shape = factor(reach), color = factor(year))) + 
    geom_point() + 
    theme_lan()

















glmm_nocov <- glmer(count ~ offset(log(area)) + flow + (flow | reach) + 
                        (1 | reach), glmerControl(calc.derivs = FALSE), 
                   data = yr_reach, family = poisson)


plot(resid(glmm_nocov) ~ yr_reach$flow)

mod_b1 <- bootMer(glmm_nocov, function(x){as.numeric(fixef(x)[2])}, nsim = 1000, 
                  type = 'parametric')
b1_ci <- as.numeric(quantile(mod_b1$t, probs = c(0.025, 0.975)))

summary(glmm_nocov)
coef(glmm_nocov)




# glmm_full <- glmer(count ~ offset(log(area)) + flow + (1 + flow | reach), 
#                    data = yr_reach, family = poisson)
# summary(glmm_full)
# full_b1 <- bootMer(glmm_full, function(x){fixef(x)[2]}, nsim = 1000)
# pchisq(2*(as.numeric(logLik(glmm_full)) - as.numeric(logLik(glmm_nocov))), df=1, 
#        lower.tail = F)
# 
# quantile(full_b1$t, probs = c(0.025, 0.975))


# glmm_int <- glmer(count ~ offset(log(area)) + flow + (1 | reach), 
#                   data = yr_reach, family = poisson)
# glmm_slp <- glmer(count ~ offset(log(area)) + flow + (0 + flow | reach), 
#                   data = yr_reach, family = poisson)
# 
# # Test the significance of the intercept random effect
# pchisq(2*(as.numeric(logLik(glmm_full)) - as.numeric(logLik(glmm_slp))), 
#                   df=2, lower.tail = F)
# [1] 8.734735e-319
# # Test the significance of the slope random effect
# pchisq(2*(as.numeric(logLik(glmm_full)) - as.numeric(logLik(glmm_int))), 
#        df=2, lower.tail = F)
# [1] 0





# Test statistics should be number of years...










