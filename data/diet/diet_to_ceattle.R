# Script for updating the hake intraspecies predation data for use in CEATTLE

library(tidyverse)
library(viridis)
library(reshape2)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Read in full aged dataset
aged_dataset <- read.csv("data/diet/CCTD_FEAT_combined.csv")

# Create overall intraspecies predation dataset -------------------------------
# Find the hake proportion for each predator
aged_wt <- aged_dataset %>%
  group_by(Predator_ID) %>%
  mutate(stomach_wt = sum(prey_wt, na.rm = TRUE)) %>%
  mutate(hake_prey_prop = if_else(prey_name == "Pacific Hake", prey_wt / stomach_wt, 0)) %>%
  select(Predator_ID, year, predator_age, prey_name, prey_age, hake_prey_prop) %>%
  distinct() # Remove duplicate rows - same pred ID, multiple hake prey 

# Total number of stomachs
stomach_n <- aged_wt %>%
  group_by(predator_age) %>%
  summarize(sample_size = n())

# Calculate average as sum of proportions per pred/prey age combo / number of stomachs per predator age
hake_prop <- aged_wt %>%
  group_by(predator_age, prey_age) %>%
  summarize(sum_prop = sum(hake_prey_prop)) %>%
  filter(!is.na(prey_age)) %>%
  left_join(stomach_n) %>%
  mutate(wt_prop = sum_prop / sample_size)

# Plot diet data
diet_plot <- ggplot(hake_prop, aes(x=as.factor(predator_age), y=wt_prop, fill=as.factor(prey_age))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(limits = factor(1:15)) +  # add in missing predator ages
  xlab("predator hake age") + ylab("diet proportion by weight") +
  labs(fill = "prey hake age")
diet_plot

ggsave(filename = "plots/diet/cannibalism_overall.png", diet_plot, 
       bg = "transparent", width=200, height=120, units="mm", dpi=300)


### Get data ready to be added directly to CEATTLE ----------------------------
# Create empty dataframe in shape for CEATTLE
all_ages <- data.frame(predator_age = rep(1:15, each = 15),
                       prey_age = rep(1:15, 15),
                       sample_size = rep(NA, 225),
                       wt_prop = rep(NA, 225))

# Merge with data
to_ceattle <- hake_prop %>%
  select(predator_age, prey_age, sample_size, wt_prop) %>%
  full_join(all_ages) %>%
  arrange(predator_age, prey_age) %>%
  distinct(predator_age, prey_age, .keep_all = TRUE) %>%
  fill(sample_size)  # sometimes need to restart R for this line to work
to_ceattle[is.na(to_ceattle)] <- 0

intrasp_ceattle <- cbind(Pred = rep(1, nrow(to_ceattle)),
                         Prey = rep(1, nrow(to_ceattle)),
                         Pred_sex = rep(0, nrow(to_ceattle)),
                         Prey_sex = rep(0, nrow(to_ceattle)),
                         Pred_age = to_ceattle$predator_age,
                         Prey_age = to_ceattle$prey_age,
                         Year = rep(0, nrow(to_ceattle)),
                         Sample_size = to_ceattle$sample_size,
                         Stomach_proportion_by_weight = to_ceattle$wt_prop)

write.csv(intrasp_ceattle, "data/diet/diet_for_CEATTLE_original.csv", row.names = FALSE)


### See if it's worth doing time-varying (yearly) predation -------------------
# same process for calculating the average as above, just disaggregated by year
stomach_n_yearly <- aged_wt %>%
  group_by(year, predator_age) %>%
  summarize(sample_size = n())
