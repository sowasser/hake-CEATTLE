# Script for combining the SWFSC/CCTD data and the NWFSC data from Alicia 
# together for Pacific Hake

library(dplyr)

CCTD_pred <- read.csv("data/diet/Full dataset/hake_aged_pred.csv")
CCTD_prey <- read.csv("data/diet/Full dataset/hake_aged_prey.csv")
# Combine predators and prey together
CCTD_all <- merge(CCTD_pred, CCTD_prey, by = "Predator_ID") %>%
  select(Month, Year, Latitude, Longitude, pred_ages, Prey_Com_Name, 
         Prey_Weight_g, prey_ages)

# Add column for source of data
CCTD_all <- cbind(CCTD_all, source = rep("CCTD", length(CCTD_all[, 1])))

###### REMOVE YEARS COVERED BY NWFSC ######


NWFSC_all <- read.csv("data/diet/all_FEAT_diet.csv")
NWFSC_all$tow_timestamp <- as.Date(NWFSC_all$tow_timestamp)
NWFSC_all$month <- lubridate::month(lubridate::ymd(NWFSC_all$tow_timestamp))
NWFSC_all$year <- lubridate::year(lubridate::ymd(NWFSC_all$tow_timestamp))

NWFSC_all <- NWFSC_all %>%
  select(month, year, tow_latitude, tow_longitude, predator_age, 
         prey_category_long, content_wt_g)

NWFSC_all <- cbind(NWFSC_all, prey_ages = rep(1, length(NWFSC_all[, 1])),
                   source = rep("NWFSC", length(NWFSC_all[, 1])))

labels <- c("month", "year", "latitude", "longitude", "predator_age", 
            "prey_name", "prey_wt", "prey_age", "source")

colnames(CCTD_all) <- labels
colnames(NWFSC_all) <- labels

all_data <- rbind(CCTD_all, NWFSC_all)

###### CHANGE PREY NAMES TO PACIFIC HAKE / OTHER #######
###### REMOVE AGES FOR PREY THAT ISN'T HAKE ######
