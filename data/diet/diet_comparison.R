# Script for comparison of hake diet datasets
library(tidyverse)

### Read in & investigate SWFSC dataset ---------------------------------------
# Read in collection info for dataset 3 - NWFSC data from Alicia 
SWFSC_collect3 <- read.csv("data/diet/Full dataset/v4/collection_information_v4.csv") %>%
  filter(Data_Set == 3)

SWFSC_predator <- read.csv("data/diet/Full dataset/v4/predator_information_v4.csv")

# Save only predators where collection IDs were found in dataset 3 
SWFSC_pred3 <- subset(SWFSC_predator, Collection_ID %in% SWFSC_collect3$Collection_ID) 

# Check that all predators are hake (they are!)
unique(SWFSC_pred3$Predator_Com_Name)

# Read in and save only prey that match the predator IDs in dataset 3
SWFSC_prey3 <- read.csv("data/diet/Full dataset/v4/prey_composition_v4.csv") %>%
  filter(Predator_ID %in% SWFSC_pred3$Predator_ID)

# Check for instances of cannibalism
SWFSC_cannibals <- SWFSC_prey3 %>%
  filter(Prey_Com_Name == "Pacific Hake")  # None!!

unique(SWFSC_prey3$Prey_Com_Name)  # no gadoids, the generalized category from NWFSC

SWFSC_gadoids <- SWFSC_prey3 %>%
  filter(Prey_Com_Name == "bony fish")  # closest to hake?


# Look a little at predator structure
unique(SWFSC_collect3$Year)  # All years represented


### Read in and investigate NWFSC dataset -------------------------------------
NWFSC_full <- read.csv("data/diet/old/full_stomachs2.csv")


# Compare two dataframes ------------------------------------------------------
SWFSC <- as.data.frame(cbind(latitude = SWFSC_collect3$Latitude,
                             longitude = SWFSC_collect3$Longitude))

NWFSC <- as.data.frame(cbind(latitude = unique(NWFSC_full$td_latitude),
                             longitude = unique(NWFSC_full$td_longitude)))

# Find if values in NWFSC cannibalism data are in the full SWFSC dataset
matched <- do.call(paste0, NWFSC) %in% do.call(paste0, SWFSC)
length(which(matched))  # 10 / 16 are....
