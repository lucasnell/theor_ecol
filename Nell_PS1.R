
# Install required packages if they're not already installed, then load
for (f in c('dplyr', 'ggplot2', 'tidyr', 'broom', 'lme4')) {
    if (!f %in% rownames(installed.packages())) {
        install.packages(f, dependencies = TRUE)
    }
    library(f, character.only = TRUE)
}; rm(f)



# =========================
# Assignment:
# =========================
# 1. Upload the grouse data and take a look at it. Zoo540_2016_PS1.R has some simple 
# plots. In addition to looking at the entire data set, compare the three different 
# PERIODS and look at the three common grouse species separately (RUGR, WITU, and STGR).
# 2. Perform one or more analyses to ask whether the effect of windspeed (either the 
# windspeed measured at individual stations or the mean windspeed within routes) 
# affects the number of RUGR during PERIOD = 3.
# 3. Prepare to demonstrate your analyses in class if you have one or more that you are 
# happy with.


# =========================
# Reading data
# =========================

dat <- read.csv(file = "grouse_data_7Sep16.csv", header = TRUE) %>%
    as.tbl %>%
    mutate(DATE = as.Date(DATE, format = "%d-%b")) %>%
    mutate_each(funs(as.factor), STATION, PERIOD)

# dplyr::as.tbl makes printing data frames less of a hazard
dat



# =========================
# Visualizing entire dataset
# =========================

# "Tidy" data for easier plotting & analyzing
grouse_df <- dat %>%
    gather(species, detected, RUGR, WITU, STGR) %>%
    mutate(
        species = factor(species, levels = c('STGR', 'RUGR', 'WITU'), 
                         labels = c('Sharp-tailed Grouse', 'Ruffed Grouse', 
                                    'Wild Turkey'))
    )

# ggplot objects relevant to all plots of spatial distributions of detections
base_tile <- {
    # Calculating proportion of visits at each location that detected grouse on the fly
    grouse_df %>%
        group_by(X_NAD83, Y_NAD83, species, PERIOD) %>%
        summarize(detected = mean(detected)) %>%
        ungroup %>%
        select(species, PERIOD, X_NAD83, Y_NAD83, detected)
    } %>%
    # Now for the plot itself
    ggplot(aes(x = X_NAD83, y = Y_NAD83, fill = detected)) +
    geom_tile(aes(width = 1e3, height = 1e3)) +
    # Asthetics...
    theme_bw() + 
    scale_fill_gradient(low = "lightblue", high = "black", guide = "colorbar") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) + 
    scale_x_continuous("Longitude") + 
    scale_y_continuous("Latitude")

# Separate by species and period
base_tile + facet_grid(PERIOD ~ species)




# =========================
# Detection ~ wind speed for ruffed grouse in period 3
# =========================

ruf3 <- grouse_df %>%
    filter(species == 'Ruffed Grouse', PERIOD == 3)

# --------
# GLM
# --------
ruf3_glm <- glm(detected ~ WIND_SPEED, binomial('logit'), data = ruf3)

summary(ruf3_glm)
AIC(ruf3_glm)

# Diagnostic plots
par(mfrow=c(2,2))
plot(ruf3_glm)
par(mfrow=c(1,1))

# Variability in model predictions
set.seed(999)
boot_ruf3_glm <- ruf3 %>% 
    bootstrap(100) %>%
    do(augment(
        glm(detected ~ WIND_SPEED, binomial('logit'), data = .),
        type.predict = 'response'
    ))
# Plotting it...
ggplot(boot_ruf3_glm, aes(x = WIND_SPEED, y = detected)) + 
    geom_point() +
    geom_line(aes(y = .fitted, group = replicate), alpha = 0.15, color = 'dodgerblue') + 
    stat_smooth(method = 'glm', method.args = list(family = binomial('logit')), 
                color = 'black') +
    # Asthetics
    theme_bw() +
    xlab(expression("Wind speed (km " * hr^{-1} * ")")) + 
    ylab("Detection of ruffed grouse")



# --------
# GLMM
# --------

ruf3_glmm <- glmer(detected ~ WIND_SPEED + (1 | ROUTE/STATION), 
                   data = ruf3, family = binomial('logit'))

# Various formats of output
summary(ruf3_glmm)
AIC(ruf3_glmm)
tidy(ruf3_glmm)

# Diagnostic plots??
plot(ruf3_glmm)



# --------
# GLMM w/ spatial autocorrelation
# --------
# library('spaMM')
# ruf3_glmm_sa <- corrHLfit(detected ~ WIND_SPEED + (1 | ROUTE/STATION) + 
#                               Matern(1 | X_NAD83 + Y_NAD83), 
#                           data = ruf3, family = binomial())
# 




# # Moving this file and CSV to Box folder...
# system(
#     paste("cd", getwd(),
#           "&& cp grouse_data_7Sep16.csv Nell_ps1.R",
#           "~/'Box Sync/ZooEnt_540_2016/Homework Folders/L_Nell/ps1/'")
# )

