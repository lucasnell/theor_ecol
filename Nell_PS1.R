
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

# dplyr::as.tbl makes printing data frames easier (many cols is still a problem)
dat



# =========================
# Visualizing entire dataset
# =========================

# "Tidy" data for easier plotting
grouse_df <- dat %>%
    gather(species, detected, RUGR, WITU, STGR) %>%
    mutate(
        species = factor(species, levels = c('STGR', 'RUGR', 'WITU'), 
                         labels = c('Sharp-tailed Grouse', 'Ruffed Grouse', 
                                    'Wild Turkey')))

# ggplot objects relevant to all plots of spatial distributions of detections
base_plot <- grouse_df %>%
    ggplot(aes(x = X_NAD83, y = Y_NAD83, color = as.factor(detected))) +
    geom_point(shape = 1, size = 2, alpha = 0.5) + 
    # Asthetics...
    scale_color_manual(values = c('tan1', 'black')) +
    theme_bw() + 
    theme(legend.position = 'none', axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    scale_x_continuous("Longitude") + 
    scale_y_continuous("Latitude")

# Separate by species only
base_plot + facet_grid(~ species)

# Separate by species and period
base_plot + facet_grid(PERIOD ~ species)





# =========================
# Detection ~ wind speed for ruffed grouse in period 3
# =========================

ruf3 <- grouse_df %>%
    filter(species == 'Ruffed Grouse', PERIOD == 3)

# --------
# GLM
# --------
ruf3_glm <- ruf3 %>%
    glm(detected ~ WIND_SPEED, binomial('logit'), data = .)

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

ggplot(boot_ruf3_glm, aes(x = WIND_SPEED, y = detected)) + 
    geom_point() +
    geom_line(aes(y = .fitted, group = replicate), alpha = 0.2, color = 'dodgerblue') + 
    # Asthetics
    theme_bw() +
    xlab(expression("Wind speed (km " * hr^{-1} * ")")) + 
    ylab("Detection of ruffed grouse")


# --------
# GLMM
# --------

ruf3_glmm <- ruf3 %>%
    glmer(detected ~ WIND_SPEED + (1 | ROUTE/STATION), 
          data = ., family = binomial('logit'))

summary(ruf3_glmm)
AIC(ruf3_glmm)


# Diagnostic plots??
plot(ruf3_glmm)

