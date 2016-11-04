
# Install required packages if they're not already installed, then load
for (f in c('magrittr', 'dplyr', 'readr', 'tidyr', 'ggplot2', 'parallel', 'lme4',
            'grid', 'nlme', 'car')) {
    if (!f %in% rownames(installed.packages())) { install.packages(f) }
    library(f, character.only = TRUE)
}; rm(f)


log_min <- function(x) {
    nz <- min(x[x > 0])
    log(x + nz)
}

z_trans <- function(x, trans = NULL) {
    if (!is.null(trans)) {
        new_x <- trans(x)
    } else {
        new_x <- x
    }
    (new_x - mean(new_x)) / sd(new_x)
}



phy_data <- read_csv(paste0("./box/ZooEnt_540_2016/Data set analyses/Phytobenthos data ",
                            "31Oct16/phytobenthos_data.csv"), 
                     col_types = 'ccDiccccdddddddd', 
                     locale = locale(date_format="%m/%d/%Y")) %>% 
    mutate_at(vars(light_dark, rack, n_fact, p_fact), funs(as.factor))

#####################
# Column descriptions
#####################

# exp: denotes 'lake' or 'lab' experiment (it might be most useful to just ignore the lab
#       experiment for the purposes of class)
# id: mesocosm ID. Lake experiment is 1-72 (but 70 does not exist). Lab is L1-L54.
# samp_date: date sample was taken
# samp_time: measurements were made on mesocosms more than once during the experiment. 
#       The different time points are denoted as 0 for 't0' (only done for the chl and 
#       om measurements on the lake mesocosms), 1 for 't1', 2 for 't2'
# light_dark: Attempted manipulation of light levels (these techniques differed in the 
#       lake vs lab experiment). Lake mesocosms were placed at different depths with 
#       the goal of achieving 'light' and 'dark' conditions. In the lab 'light' and 
#       'dark' was done with shading cloth.
# rack: mesocosms were grouped into different racks for experimental deployment (or 
#       coolers for the lab experiment). Each rack contained a complete block of all 9
#       treatment combinations and 3 additional combinations that were haphazardly 
#       assigned.
# n_fact: Nitrogen treatment as a factor. Control: no NH4Cl added; Medium: 0.1 M NH4Cl 
#       added; High: 0.5 M NH4Cl added
# p_fact: Phosphorus treatment as a factor. Control: no KH2PO4 added; Medium: 0.0065 M 
#       KH2PO4 added; High: 0.033 M KH2PO4 added
# n_conc: Nitrogen -as NH4Cl- concentration (M); controls have ambient sediment 
#       concentration of NH4
# p_conc: Phosphorus -as KH2PO4- concentration (M); controls have ambient sediment 
#       concentration of PO4 
# om: percent organic matter in sediment. Calculated as: (dry weight (mg) - combusted 
#       weight (mg)) / dry weight (mg)
# chl: sediment chlorophyll a concentration (ug/l)
# chl_om: chlorophyll a concentration relative to organic matter content
# nep: Net Ecosystem Production (change in water percent O2 saturation during light 
#       incubation)
# resp: Ecosystem Respiration (change in water percent O2 saturation during dark 
#       incubation)
# gpp: Gross Primary Production (NEP+Resp; units: delta O2 saturation)

# There are occasional NA's throughout the data set. 
# Sometimes this was because something was not measured at that sample time (nep, resp, 
#   & gpp during time point '0')
# Other times a sample for a specific measurement was lost, removed from data due to lab 
#   errors, etc.



# Renaming factor levels for N and P treatments
phy_data <- phy_data %>%
    mutate(
        n_fact = recode_factor(n_fact, con = 'control N', 
                               med = 'medium N', high = 'high N'),
        p_fact = recode_factor(p_fact, con = 'control P', 
                               med = 'medium P', high = 'high P')
    )

lake_data <- phy_data %>% 
    filter(exp == 'lake', samp_time > 0) %>% 
    select(-exp, -om, -chl, -chl_om, -resp, -gpp)




lake_data %>% 
    group_by(id) %>% 
    summarize(d_nep = nep[samp_time == 2] - nep[samp_time == 1],
              # r_nep = nep[samp_time == 2] / nep[samp_time == 1],
              p_fact = p_fact[1],
              n_fact = n_fact[1],
              light_dark = light_dark[1],
              rack = rack[1]) %>% 
    ggplot(aes(light_dark, d_nep)) + 
    geom_point(position = position_jitter(width = 0.3, height = 0), alpha = 1.0) +
    stat_summary(geom = 'errorbar', fun.data = 'mean_cl_boot',
                 width = 0.2, color = 'red', size = 0.75) +
    facet_grid(n_fact ~ p_fact) +
    theme_bw() +
    ylab(expression(Delta * 'NEP')) +
    xlab('')


z <- lmer(log(nep) ~ samp_time * light_dark * n_fact * p_fact + (1 | id) + (1 | rack), 
          data = lake_data, REML = FALSE)
summary(z)

Anova(z, type = 2)





p2 <- lake_data %>% 
    group_by(id) %>% 
    summarize(d_nep = nep[samp_time == 2] - nep[samp_time == 1],
              p_fact = p_fact[1],
              n_fact = n_fact[1],
              light_dark = light_dark[1],
              rack = rack[1]) %>% 
    ggplot(aes(n_fact, d_nep)) + 
    geom_point(position = position_jitter(width = 0.3, height = 0), alpha = 0.3) +
    stat_summary(geom = 'errorbar', fun.data = 'mean_cl_boot',
                 width = 0.2, color = 'red', size = 0.75) +
    facet_grid( ~ light_dark) +
    theme_bw() +
    ylab(expression(Delta * 'NEP')) +
    xlab('N treatment')


grid.newpage(); grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = 'last'))



# p1 <- phy_data %>% 
#     filter(exp == 'lake') %>% 
#     ggplot(aes(light_dark, nep)) +
#     geom_point(position = position_jitter(width = 0.3, height = 0), alpha = 0.3) +
#     stat_summary(geom = 'errorbar', fun.data = 'mean_cl_boot',
#                  width = 0.2, color = 'red', size = 0.75) +
#     facet_grid(p_fact ~ n_fact) +
#     theme_bw()
# 
# p2 <- phy_data %>% 
#     filter(exp == 'lake') %>% 
#     ggplot(aes(p_conc, nep, color = light_dark)) +
#     geom_point(alpha = 0.3) +
#     geom_smooth(method = 'lm', se = FALSE) +
#     theme_bw()
# 
# p3 <- phy_data %>% 
#     filter(exp == 'lake') %>% 
#     ggplot(aes(n_conc, nep, color = light_dark)) +
#     geom_point(alpha = 0.3) +
#     geom_smooth(method = 'lm', se = FALSE) +
#     theme_bw()
# 
# 
# p1
# grid.newpage(); grid.draw(cbind(ggplotGrob(p2 + theme(legend.position = 'none')), 
#                                 ggplotGrob(p3), size = 'last'))
# 
# phy_data$id %>% unique %>% length
