# Script for subsetting the whole diet database for predator hake

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(forcats)

path <- "data/diet/Full dataset/v4/"

# Read in and subset datasets -------------------------------------------------
# Filter predator dataset for just hake and collections from 1980 onwards
predator_hake <- read.csv(paste0(path, "predator_information_v4.csv")) %>%
  filter(Predator_Com_Name == "Pacific Hake") %>%
  filter(Collection_ID > 58)

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


# Plot general trends in data -------------------------------------------------
pred_fl <- ggplot(predator_hake, aes(x = FL_cm)) +
  geom_histogram() +
  xlab("fork length (cm)") + ylab(" ") +
  theme_sleek()
pred_fl

prey_sp <- ggplot(prey_of_hake_comp, aes(y = fct_infreq(Prey_Com_Name))) +
  geom_bar(position = "dodge") +
  theme_sleek()
prey_sp

prey_length <- ggplot(hake_prey, aes(x = (Prey_Length1/10))) +
  geom_histogram() +
  xlab("length (cm)") + ylab(" ") +
  theme_sleek()
prey_length

