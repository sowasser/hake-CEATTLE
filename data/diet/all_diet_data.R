# Script for subsetting the whole diet database for predator hake

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(forcats)

path <- "data/diet/Full dataset/v4/"

### Read in and subset datasets -----------------------------------------------
# Filter predator dataset for just hake and collections from 1980 onwards
predator_hake <- read.csv(paste0(path, "predator_information_v4.csv")) %>%
  filter(Predator_Com_Name == "Pacific Hake") %>%
  filter(Collection_ID > 58) %>%
  filter(!is.na(FL_cm))  # Remove observations without fork lengths

# Filter collection info by hake predator collection IDs
collection_hake <- read.csv(paste0(path, "collection_information_v4.csv")) %>%
  filter(Collection_ID %in% predator_hake$Collection_ID)

# Filter prey composition dataset by hake predator IDs
prey_of_hake_comp <- read.csv(paste0(path, "prey_composition_v4.csv")) %>%
  filter(Predator_ID %in% predator_hake$Predator_ID)

# Filter prey size dataset by hake-predated-prey IDs
prey_of_hake_size <-  read.csv(paste0(path, "prey_size_v4.csv")) %>%
  filter(Prey_Comp_ID %in% prey_of_hake_comp$Prey_Comp_ID)

# Subset prey dataset by prey that are hake
hake_prey <- prey_of_hake_size %>%
  filter(Prey_Com_Name == "Pacific Hake")


### Combine then refine datasets ----------------------------------------------
# Combine collection info with predator info
all_pred <- merge(collection_hake, predator_hake, all.y = TRUE)[, c("Month", 
                                                                    "Year", 
                                                                    "Latitude", 
                                                                    "Longitude",
                                                                    "Predator_ID", 
                                                                    "FL_cm", 
                                                                    "Predator_Weight_kg")]

# Combine prey-of-hake datasets together
all_prey <- merge(prey_of_hake_comp, prey_of_hake_size, all.y = TRUE)[, c("Prey_Com_Name", 
                                                                          "Predator_ID", 
                                                                          "Prey_Weight_g",
                                                                          "Prey_Maturity",
                                                                          "Prey_Length1")]


# Subset cannibalized prey
hake_prey <- all_prey %>%
  filter(Prey_Com_Name == "Pacific Hake")

# Combine pred & prey togehter - this is just full hake stomachs.
full_stomachs <- merge(all_pred, all_prey, all.y = TRUE)

# Full cannibalism instances
hake_hake <- merge(all_pred, hake_prey, all.y = TRUE)


### Plot general trends in data -----------------------------------------------
pred_fl <- ggplot(predator_hake, aes(x = FL_cm)) +
  geom_histogram() +
  xlab("fork length (cm)") + ylab(" ") +
  theme_sleek()
pred_fl

prey_sp <- ggplot(prey_of_hake_comp, aes(y = fct_infreq(Prey_Com_Name))) +
  geom_bar(position = "dodge") +
  xlab(" ") + ylab("prey common name") +
  theme_sleek()
prey_sp

prey_length <- ggplot(hake_prey, aes(x = (Prey_Length1/10))) +
  geom_histogram() +
  xlab("length (cm)") + ylab(" ") +
  theme_sleek()
prey_length


### Combine predator & prey datasets & write new .csvs ------------------------
write.csv(all_pred, "data/diet/full_hake_pred.csv")
write.csv(all_prey, "data/diet/full_prey.csv")

all_hake <- merge(all_pred, all_prey, all = TRUE)
write.csv(all_hake, "data/diet/full_hake_diet.csv")
