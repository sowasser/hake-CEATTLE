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
  filter(Collection_ID > 58) %>%  # Years before 1980
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


# Combine pred & prey togehter - this is just full hake stomachs.
full_stomachs <- merge(all_pred, all_prey, all.y = TRUE)

# Cannibalism instances
hake_hake <- full_stomachs %>%
  filter(Prey_Com_Name == "Pacific Hake")


### Plot general trends in data -----------------------------------------------
# Occurrences of predation over 100 n
high_wt <- prey_of_hake_comp %>% 
  group_by(Prey_Com_Name) %>% 
  summarize(highest = sum(Prey_Weight_g)) %>%
  slice_max(n = 10, order_by = highest)

high_n <- prey_of_hake_comp %>% 
  group_by(Prey_Com_Name) %>% 
  summarize(highest = n()) %>%
  slice_max(n = 10, order_by = highest)

highest <- rbind(cbind(high_wt, variable = rep("weight (top 10)", 10)), 
                 cbind(high_n, variable = rep("occurrence (top 10)", 10)))

prey_sp <- ggplot(highest, aes(x = Prey_Com_Name, y = highest,
                               fill = ifelse(Prey_Com_Name == "Pacific Hake", "highlighted", "normal"))) +
  geom_bar(position = "dodge", stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
  theme_sleek() +
  xlab(" ") + ylab(" ") +
  facet_wrap(~ variable, scales = "free")
prey_sp

ggsave(filename = "plots/diet/hake_prey_species.png", prey_sp, 
       width=300, height=150, units="mm", dpi=300)

# Lengths
pred_fl <- ggplot(predator_hake, aes(x = FL_cm)) +
  geom_histogram() +
  xlab("fork length (cm)") + ylab(" ") +
  theme_sleek()
pred_fl

prey_length <- ggplot(hake_hake, aes(x = (Prey_Length1/10))) +
  geom_histogram() +
  xlab("length (cm)") + ylab(" ") +
  theme_sleek()
prey_length


# Timing
timing_all <- rbind(cbind(all_pred[, c("Month", "Year")], type = rep("all predator hake", length(all_pred$Month))),
                    cbind(hake_hake[, c("Month", "Year")], type = rep("cannibalistic hake", length(hake_hake$Month)))) %>%
  group_by(Year, Month, type) %>%
  summarize(n = n()) %>%
  filter(!is.na(Year)) 

##### NB: This plot is not quite correct! Instances of cannibalism are counted twice #####
timing <- ggplot(timing_all, aes(x = as.factor(Month), y = n, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  facet_wrap(~ Year)
timing
  

### Write predator & prey datasets to .csvs -----------------------------------
write.csv(all_pred, "data/diet/Full dataset/full_hake_pred.csv")
write.csv(all_prey, "data/diet/Full dataset/full_prey.csv")
