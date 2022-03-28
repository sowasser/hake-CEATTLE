# Script for manipulating the ROMS data - two different products from 1980 - 
# 2010 and 2011-2020 - to work as both a timeseries for the CEATTLE model and
# for parameterizing the hake bioenergetics (found in hake_temperature.R)

library(dplyr)

# Read in ROMS datasets & keep the columns that actually include data
ROMS_8010 <- read.csv("data/temperature/ROMS_temp_hake_south_1980_2010.csv")[, c(1:4)]
ROMS_1120 <- read.csv("data/temperature/ROMS_temp_hake_south_2011_2020.csv")[, c(1:4)]

ROMS_all <- rbind(ROMS_8010, ROMS_1120)

# Find yearly means and export 
ROMS_all_mean <- ROMS_all %>% 
  group_by(Year) %>%
  summarise(mean_temp = mean(H2)) 

write.csv(ROMS_all_mean, "data/temperature/ROMS_mean.csv", row.names = FALSE)

# Isolate summer months, find yearly means, and export
ROMS_summer <- ROMS_all %>% 
  filter(between(Month, 6, 9)) %>%
  group_by(Year) %>%
  summarise(mean_temp = mean(H2))

write.csv(ROMS_summer, "data/temperature/ROMS_summer_mean.csv", row.names = FALSE)
