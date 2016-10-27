library(tidyverse)
library(dataRetrieval)

if (grepl('theor_ecol$', getwd())){
    setwd('./Pine_dat')
} else if (!grepl('Pine_dat$', getwd())) {
    stop('Unknown directory. It should be in one ending in "theor_ecol" or "Pine_dat".')
}

# ------------
# Angostura and San Acacia flow
# ------------
# Use these parameters to download the station data
siteNumbers <- c('08330000', '08358400')
# 08330000 = RIO GRANDE AT ALBUQUERQUE, NM
# 08358400 = RIO GRANDE FLOODWAY AT SAN MARCIAL, NM
# For Angostura and San Acacia reaches, resp.
parameterCd <- "00060"        # discharge cfs
dailyMean <- "00003"        # daily mean value
startDate <- "1993-01-01"
endDate <- "2015-12-31"

# This will show what is available for the station
# whatNWISdata(siteNumbers, service="dv", statCd=dailyMean)


# Pull in the data for the parameters above
flowD <- readNWISdv(siteNumbers,parameterCd, startDate, endDate, statCd=dailyMean)
#  ... also can request  statCd=statCd
flowD <- renameNWISColumns(flowD)
flowD <- as.tbl(flowD) %>% mutate(
    reach = ifelse(site_no == '08330000', 'Angostura', 
                   ifelse(site_no == '08330875', 'Isleta', 
                          'San Acacia')),
    year = Date %>% format(., '%Y') %>% as.integer,
    month = Date %>% format(., '%m') %>% as.integer) %>% 
    filter(month %in% c(5,6))

# ------------
# Isleta flow
# ------------
# 08331510 = RIO GRANDE AT STATE HWY 346 NEAR BOSQUE, NM
# 08332010 = RIO GRANDE FLOODWAY NEAR BERNARDO, NM
# For Isleta reach, 2006-present and 1964-2005, respectively
# (The second one has no data for 2006-2011)
whatNWISdata('08332010', service="dv", statCd=dailyMean)

old_f <- readNWISdv('08331510', parameterCd, '2006-01-01', endDate, statCd=dailyMean)
old_f <- renameNWISColumns(old_f)

new_f <- readNWISdv('08332010', parameterCd, startDate, '2005-12-31', statCd=dailyMean)
new_f <- renameNWISColumns(new_f)

flowD_Isl <- list(old_f, new_f) %>% bind_rows %>% as.tbl %>% 
    mutate(
    reach = 'Isleta',
    year = Date %>% format(., '%Y') %>% as.integer,
    month = Date %>% format(., '%m') %>% as.integer) %>% 
    filter(month %in% c(5,6))





# left off: assemble these below





duration_flowD <- bind_rows(flowD, flowD_Isl) %>% 
    group_by(year, reach) %>% 
    summarize(
        mean = mean(Flow) %>% round(., 2),
        sd = sd(Flow) %>% round(., 2),
        amp = diff(range(Flow)) %>% as.integer,
        over_1000 = length(Flow[Flow > 1000]),
        over_2000 = length(Flow[Flow > 2000]),
        over_3000 = length(Flow[Flow > 3000]),
        over_4000 = length(Flow[Flow > 4000]),
        over_5000 = length(Flow[Flow > 5000])
    ) %>% 
    ungroup




write_csv(duration_flowD, 'all_flow_data.csv')





