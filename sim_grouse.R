
# Install required packages if they're not already installed, then load
for (f in c('dplyr', 'ggplot2', 'grid', 'parallel', 'tidyr')) {
    if (!f %in% rownames(installed.packages())) {
        install.packages(f, dependencies = TRUE)
    }
    library(f, character.only = TRUE)
}; rm(f)


# 5-8 stations, each of 178 routes
# obs period = 3 min


# =========================
# Reading data
# =========================

dat <- read.csv(file = "grouse_data_7Sep16.csv", header = TRUE) %>%
    as.tbl %>%
    mutate(DATE = as.Date(DATE, format = "%d-%b")) %>%
    mutate_each(funs(as.factor), STATION, PERIOD)


# "Tidy" data for easier plotting & analyzing
grouse_df <- dat %>%
    gather(species, detected, RUGR, WITU, STGR) %>%
    mutate(
        species = factor(species, levels = c('STGR', 'RUGR', 'WITU'), 
                         labels = c('Sharp-tailed Grouse', 'Ruffed Grouse', 
                                    'Wild Turkey'))
    )

asinh(grouse_df$WIND_SPEED) %>% hist(.)

# model for variation among routes
# another for within routes



# More variation by routes than predicted by if they had the same among routes

