# Script for combining the SWFSC/CCTD data and the FEAT data from Alicia 
# together for Pacific Hake

library(dplyr)

### Read in & update CCTD data from SWFSC -------------------------------------
CCTD_pred <- read.csv("data/diet/Full dataset/hake_aged_pred.csv")
CCTD_prey <- read.csv("data/diet/Full dataset/hake_aged_prey.csv")

# Combine predators and prey together, generalize prey names to "other"
CCTD_all <- merge(CCTD_pred, CCTD_prey, by = "Predator_ID") %>%
  mutate(prey_name = ifelse(Prey_Com_Name == "Pacific Hake", "Pacific Hake", "other")) %>%
  select(Month, Year, Latitude, Longitude, pred_ages, prey_name, 
         Prey_Weight_g, Prey_Length1, prey_ages) %>%
  filter(Year <= 2004)  # remove years covered by FEAT 

# Add column for source of data
CCTD_all <- cbind(CCTD_all, source = rep("CCTD", length(CCTD_all[, 1])))


### Read in & update FEAT data from NWFSC -------------------------------------
FEAT_data <- read.csv("data/diet/all_FEAT_diet.csv")
FEAT_prey_lengths <- read.csv("data/diet/Old/hake_prey_lengths_FEAT.csv")

# Add prey lengths based on stomach ID
FEAT_all <- merge(FEAT_data, FEAT_prey_lengths, by = "stomach_uuid", all.x = TRUE)

# Expand to month and year columns
FEAT_all$tow_timestamp <- as.Date(FEAT_all$tow_timestamp)
FEAT_all$month <- lubridate::month(lubridate::ymd(FEAT_all$tow_timestamp))
FEAT_all$year <- lubridate::year(lubridate::ymd(FEAT_all$tow_timestamp))

# Generalize prey categories to Pacific Hake & "other", select columns
FEAT_all <- FEAT_all %>%
  mutate(prey_name = ifelse(prey_category_long == "Gadiformes", "Pacific Hake", "other")) %>%
  select(month, year, tow_latitude, tow_longitude, predator_age, 
         prey_name, content_wt_g, measure_value)
  

### Combine datasets ----------------------------------------------------------
FEAT_all <- cbind(FEAT_all, prey_ages = rep(1, length(FEAT_all[, 1])),
                   source = rep("FEAT", length(FEAT_all[, 1])))

labels <- c("month", "year", "latitude", "longitude", "predator_age", 
            "prey_name", "prey_wt", "prey_length", "prey_age", "source")

colnames(CCTD_all) <- labels
colnames(FEAT_all) <- labels

all_data <- rbind(CCTD_all, FEAT_all)

# Remove all prey ages for prey that isn't hake
all_data$prey_age[all_data$prey_name == "other"] <- NA

hake_data <- all_data %>%
  filter(prey_name == "Pacific Hake")
