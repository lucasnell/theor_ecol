d <- read.csv(file="./box/ZooEnt_540_2016/Data set analyses/Phytobenthos data 31Oct16/phytobenthos_data.csv", header=T)
d$samp_date <- as.Date(d$samp_date, format="%m/%d/%Y") #excel inputs dates as 31-Oct-15, b:three letter abbrev, y:2digit year

#####################
# Column descriptions
#####################

# exp: denotes 'lake' or 'lab' experiment (it might be most useful to just ignore the lab experiment for the purposes of class)
# id: mesocosm ID. Lake experiment is 1-72 (but 70 does not exist). Lab is L1-L54.
# samp_date: date sample was taken
# samp_time: measurements were made on mesocosms more than once during the experiment. The different time points are denoted as 0 for 't0' (only done for the chl and om measurements on the lake mesocosms), 1 for 't1', 2 for 't2'
# light_dark: Attempted manipulation of light levels (these techniques differed in the lake vs lab experiment). Lake mesocosms were placed at different depths with the goal of achieving 'light' and 'dark' conditions. In the lab 'light' and 'dark' was done with shading cloth.
# rack: mesocosms were grouped into different racks for experimental deployment (or coolers for the lab experiment). Each rack contained a complete block of all 9 treatment combinations and 3 additional combinations that were haphazardly assigned.
# n_fact: Nitrogen treatment as a factor. Control: no NH4Cl added; Medium: 0.1 M NH4Cl added; High: 0.5 M NH4Cl added
# p_fact: Phosphorus treatment as a factor. Control: no KH2PO4 added; Medium: 0.0065 M KH2PO4 added; High: 0.033 M KH2PO4 added
# n_conc: Nitrogen -as NH4Cl- concentration (M); controls have ambient sediment concentration of NH4
# p_conc: Phosphorus -as KH2PO4- concentration (M); controls have ambient sediment concentration of PO4 
# om: percent organic matter in sediment. Calculated as: (dry weight (mg) - combusted weight (mg)) / dry weight (mg)
# chl: sediment chlorophyll a concentration (ug/l)
# chl_om: chlorophyll a concentration relative to organic matter content
# nep: Net Ecosystem Production (change in water percent O2 saturation during light incubation)
# resp: Ecosystem Respiration (change in water percent O2 saturation during dark incubation)
# gpp: Gross Primary Production (NEP+Resp; units: delta O2 saturation)

# There are occasional NA's throughout the data set. 
# Sometimes this was because something was not measured at that sample time (nep, resp, & gpp during time point '0')
# Other times a sample for a specific measurement was lost, removed from data due to lab errors, etc.

library(lme4)
library(lmerTest)
library(nlme)
library(car)


levels(d$n_fact)
d$n_fact <- relevel(d$n_fact,3)
levels(d$n_fact)
d$n_fact <- relevel(d$n_fact,2)
levels(d$n_fact)

d$p_fact <- relevel(d$p_fact,3)
d$p_fact <- relevel(d$p_fact,2)
levels(d$p_fact)

dd <- d[d$exp == 'lake' & d$samp_time > 0,]
dd <- dd[!is.na(dd$gpp),]

par(mfrow=c(2,1))
hist(dd$nep, breaks = 30)
hist(log(dd$nep), breaks = 30)

par(mfrow=c(2,2))
for(st in 1:2) for(ld in levels(dd$light_dark)){
	ds <- dd[dd$samp_time == st & dd$light_dark == ld,]
	interaction.plot(ds$n_fact, ds$p_fact, log(ds$nep), main=paste('time =',st, ', trt =',ld), xlab = 'N', ylab = 'log NEP', ylim = c(0,5), col = as.numeric(ds$p_fact), lty = 1)
	points(log(ds$nep) ~ as.numeric(ds$n_fact), col = as.numeric(ds$p_fact))
}

z <- lmer(log(nep) ~ samp_time * light_dark * n_fact * p_fact + (1 | id) + (1 | rack), data=dd, REML = F)
Anova(z, type = 2)

# type II or type III Wald tests?
Anova(lm(log(nep) ~ n_fact * p_fact , data=dd), type = 3)
Anova(lm(log(nep) ~ n_fact * p_fact , data=dd), type = 2)
Anova(lm(log(nep) ~ n_fact + p_fact , data=dd))

# REML or ML?
z <- lmer(log(nep) ~ samp_time * light_dark * n_fact * p_fact + (1 | id) + (1 | rack), data=dd, REML = F)
Anova(z, type = 2)

z.REML <- lmer(log(nep) ~ samp_time * light_dark * n_fact * p_fact + (1 | id) + (1 | rack), data=dd, REML = T)
Anova(z.REML, type = 2)

# model selection by p-values
z <- lmer(log(nep) ~ samp_time + light_dark + n_fact + p_fact + samp_time:light_dark + samp_time:n_fact + light_dark:n_fact + samp_time:p_fact + light_dark:p_fact + n_fact:p_fact + samp_time:light_dark:n_fact + samp_time:light_dark:p_fact + samp_time:n_fact:p_fact + light_dark:n_fact:p_fact + samp_time:light_dark:n_fact:p_fact + (1 | id) + (1 | rack), data=dd, REML = F)
Anova(z, type = 2)


z_reml <- lmer(log(nep) ~ samp_time + light_dark + n_fact + p_fact + samp_time:light_dark + samp_time:n_fact + light_dark:n_fact + samp_time:p_fact + light_dark:p_fact + n_fact:p_fact + samp_time:light_dark:n_fact + samp_time:light_dark:p_fact + samp_time:n_fact:p_fact + light_dark:n_fact:p_fact + samp_time:light_dark:n_fact:p_fact + (1 | id) + (1 | rack), data=dd, REML = TRUE)
# Resources: I put three pdfs in the Phytobenthos folder that I stole from the web:

# a. "type 2 and 3 SS.pdf" gives a summary of types of sums-of-squares. The same arguments apply to types of tests for significance.

# b. "ML vs. REML.pdf" gives a comparison between ML and REML. It is a little technical, but I thought it explained things well.

# c. "Aho_et_al-2014-Ecology.pdf" gives a description of AIC and BIC used for model selection. Please read this for Monday.




# Question 1: Perform stepwise backwards regression using both ML and REML. Do the final
# models differ? At one point in the stepwise process do they become different?

z_step <- step(z, data=dd, direction="backward")

z_reml_step <- step(z_reml, data=dd, direction="backward")


# as.data.frame(fixef(z_step$model))
# as.data.frame(fixef(z_reml_step$model))



# Question 2: Analyze resp in the same way as nep. Do you get similar models? Are the 
# ML and REML models the same?

# model selection by p-values
z2 <- lmer(log(resp) ~ samp_time + light_dark + n_fact + p_fact + samp_time:light_dark + samp_time:n_fact + light_dark:n_fact + samp_time:p_fact + light_dark:p_fact + n_fact:p_fact + samp_time:light_dark:n_fact + samp_time:light_dark:p_fact + samp_time:n_fact:p_fact + light_dark:n_fact:p_fact + samp_time:light_dark:n_fact:p_fact + (1 | id) + (1 | rack), data=dd, REML = F)
# Anova(z, type = 2)

z_reml2 <- lmer(log(resp) ~ samp_time + light_dark + n_fact + p_fact + samp_time:light_dark + samp_time:n_fact + light_dark:n_fact + samp_time:p_fact + light_dark:p_fact + n_fact:p_fact + samp_time:light_dark:n_fact + samp_time:light_dark:p_fact + samp_time:n_fact:p_fact + light_dark:n_fact:p_fact + samp_time:light_dark:n_fact:p_fact + (1 | id) + (1 | rack), data=dd, REML = TRUE)

z_step2 <- step(z2, data=dd, direction="backward")

z_reml_step2 <- step(z_reml2, data=dd, direction="backward")


as.data.frame(fixef(z_step2$model))

as.data.frame(fixef(z_reml_step2$model))

# Question 3: Graph the data to try to see evidence for the patterns implied by the best-
# fitting model. The best type of figure would be to map the fitted model onto the data,
# although this isn't simple (unless dplyr gives an automatic function for this: ask 
# Lucas).




library(tidyverse)
library(broom)

ml_vars <- strsplit((z_step$model %>% formula %>% paste)[3], ' \\+ ')[[1]]
ml_vars <- ml_vars[!grepl('[:|]', ml_vars)]

reml_vars <- strsplit((z_reml_step$model %>% formula %>% paste)[3], ' \\+ ')[[1]]
reml_vars <- reml_vars[!grepl('[:|]', reml_vars)]

all_vars <- unique(c(ml_vars, reml_vars))

dd <- dd %>% as.tbl

var_list <- list()
for (v in all_vars){
    var_list[[v]] <- unique(dd[[v]])
}

nd <- expand.grid(var_list)

nd$log_nep <- predict(z_step$model, re.form = NA, newdata = nd)
# nd$reml <- predict(z_reml_step$model, re.form = NA, newdata = nd)
# nd <- nd %>% gather(model,log_resp, ml, reml) %>% as.tbl

nd
dd

# Question 4: Graph the data for RESP along with the best-fitting model. The best type 
# of figure would be to map the fitted model onto the data, although this isn't simple. 
# Here is example code for the NEP data:


dd %>% ggplot(aes(as.integer(n_fact), nep, color = p_fact)) + 
    theme_bw() +
    geom_point(shape = 1) + 
    facet_grid(light_dark ~ samp_time) +
    stat_summary(data = nd, aes(y = exp(log_nep)), geom = 'line', fun.y = mean, 
                 size = 0.75) +
    scale_y_continuous('NEP', trans = 'log', breaks = c(1, 5, 20, 50)) +
    scale_x_continuous('N treatment', expand = c(0.1, 0), breaks = seq(3), 
                       labels = levels(dd$n_fact)) +
    scale_color_brewer(type = 'qual', palette = 'Dark2')




# Question 5: Perform model selection using AIC and BIC on the RESP data. How does the 
# suite of best-fitting models compare to what you would infer from model selection 
# using p-values? How would the best-fitting models compare to the full model with 
# p-values computed using type II Wald tests?

library(MASS)
library(lmerTest)
library(stringr)

resp_lmer <- lmer(log(resp) ~ samp_time + light_dark + n_fact + p_fact + samp_time:light_dark + samp_time:n_fact + light_dark:n_fact + samp_time:p_fact + light_dark:p_fact + n_fact:p_fact + samp_time:light_dark:n_fact + samp_time:light_dark:p_fact + samp_time:n_fact:p_fact + light_dark:n_fact:p_fact + samp_time:light_dark:n_fact:p_fact + (1 | id) + (1 | rack), data=dd, REML = F)
# Anova(z, type = 2)

y_var <- 'log(resp) ~'
me_vars <- '+ (1 | id) + (1 | rack)'


x_vars <- c("samp_time", 
"light_dark",
"n_fact", 
"p_fact")

# # 2-way interactions
# c("samp_time:light_dark", 
# "samp_time:n_fact",
# "light_dark:n_fact", 
# "samp_time:p_fact",
# "light_dark:p_fact", 
# "n_fact:p_fact")
# # 3-way interactions
# c("samp_time:light_dark:n_fact", 
# "samp_time:light_dark:p_fact",
# "samp_time:n_fact:p_fact", 
# "light_dark:n_fact:p_fact")
# # 4-way
# c("samp_time:light_dark:n_fact:p_fact")


one_size <- function(i, x_vars){
    var_matrix <- t(combn(x_vars, i))
    additives <- apply(var_matrix, 1, function(x) paste0(x, collapse = ' + '))
    if (i == 1) {
        return(additives)
    }
    interact_list <- lapply(2:i, 
                            function(j) {
                                apply(t(combn(x_vars, j)), 1, 
                                      function(x) paste0(x, collapse = ':'))})
    all_inters <- c(interact_list, recursive = TRUE)
    
    # Which interactions have terms present in each row of var_matrix
    interactions <- lapply(1:nrow(var_matrix), function(i){
        a <- var_matrix[i,]
        int_index <- sapply(all_inters, 
                            function(s){
                                ss <- strsplit(s, ':')[[1]]
                                return(all(ss %in% a))
                            }, USE.NAMES = FALSE)
        int_a <- all_inters[int_index]
        return(int_a)
    })
    
    inter_combos <- lapply(interactions, function(ints){
        c(sapply(1:length(ints), function(j){
            max_ints <- ifelse(j < 4, 2, 3)
            int_matrix <- as.matrix(t(combn(interactions[[1]], j)))
            par_nums <- apply(int_matrix, c(1,2), str_count, pattern = ':') + 1
            int_matrix_filt <- matrix(int_matrix[(rowSums(par_nums > max_ints) == 0), ], 
                                      ncol = ncol(int_matrix))
            apply(int_matrix_filt, 1, 
                  function(s) {
                      paste0(s, collapse = ' + ') 
                  })
        }), recursive = TRUE)
    })
    
    i_combos <- c(additives, lapply(1:length(additives),
                                    function(j) {
                                        adds <- additives[j]
                                        inters <- inter_combos[[j]]
                                        paste(adds, inters, sep = ' + ')
                                    }), recursive = TRUE)
    return(i_combos)
}

all_combos <- c(lapply(1:length(x_vars), one_size, x_vars = x_vars), recursive = TRUE)
all_combos %>% length

# # To find strings with both of a vector c('ee', 'dd')
# Reduce(`&`, lapply(c('ee', 'dd'), grepl, c('dd:ee', 'ff:dd', 'aa:ee')))

fits <- lapply(all_combos, function(c) {
    mod <- suppressMessages(lmer(as.formula(paste(y_var, c, me_vars)), data=dd, 
                                 REML = F))
    aic_bic <- c(AIC(mod), BIC(mod))
    return(aic_bic)
})

fit_df <- do.call(rbind, fits) %>% as_data_frame
colnames(fit_df) <- c('aic', 'bic')
fit_df <- fit_df %>% 
    mutate(fixed = all_combos,
           d_aic = aic - min(aic),
           d_bic = bic - min(bic))

fit_df %>% dplyr::select(d_aic, fixed) %>% arrange(d_aic)

fit_df %>% dplyr::select(d_bic, fixed) %>% arrange(d_bic)


AIC(lmer(log(resp) ~ samp_time + light_dark + n_fact + samp_time:n_fact + (1 | id) + (1 | rack), data = dd, REML = F))
AIC(lmer(log(resp) ~ samp_time + light_dark + n_fact + p_fact + samp_time:n_fact + samp_time:p_fact + (1 | id) + (1 | rack), data = dd, REML = F))


read_csv(file="./box/ZooEnt_540_2016/Data set analyses/Phytobenthos data 31Oct16/phytobenthos_data.csv") %>% 
    filter(samp_time == 0) %>% #group_by(rack) %>% 
    # summarize(chl_om = mean(chl_om,na.rm = TRUE))
    ggplot(aes(rack, chl_om, color = rack)) + 
    geom_point() + 
    stat_summary(geom = 'point', fun.y = 'mean', size = 4, shape = 1) +
    theme_bw() +
    theme(legend.position = 'none')





# # Moving this file to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp Nell_PS2_Phytobenthos.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/'")
# )