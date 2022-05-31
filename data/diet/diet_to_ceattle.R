# Script for updating the hake intraspecies predation data for use in CEATTLE

library(ggplot2)
library(ggsidekick)
library(viridis)
library(dplyr)
library(tidyr)

# Read in full aged datasets & remove unecessary index column
new_pred <- read.csv("data/diet/Full dataset/full_aged_pred.csv")[, -1]
new_prey <- read.csv("data/diet/Full dataset/full_aged_prey.csv")[, -1]


### Combine into new dataset --------------------------------------------------
aged_dataset <- merge(new_pred, new_prey, by = "Predator_ID")  # using all prey data here

# Look at instances where prey hake age = NA - all are immature
aged_dataset %>%
  filter(Prey_Com_Name == "Pacific Hake") %>%
  filter(is.na(prey_ages))

# Replace those NAs with age 1
aged_dataset$prey_ages[is.na(aged_dataset$prey_ages) & aged_dataset$Prey_Com_Name == "Pacific Hake"] <- 1

# Create new, organized dataframe for rates of cannibalism, fill in predator info
aged_subset <- aged_dataset[, c("Predator_ID", "Year", "pred_ages", "Prey_Com_Name", "prey_ages", "Prey_Weight_g")] %>%
  arrange(Predator_ID, Year) %>%
  fill(Year) %>%
  fill(pred_ages)


# Create overall intraspecies predation dataset -------------------------------
# Get number of stomachs per predator age
stomachs_n <- aged_subset %>%
  group_by(pred_ages) %>%
  summarize(sample_size = n())
  
# Total stomach weight per predator age
total_wt <- aged_subset %>%
  group_by(pred_ages) %>%
  summarize(total_wt = sum(Prey_Weight_g, na.rm = TRUE))

# Combine summarized datasets and calculate stomach weight proportions (overall)
intrasp <- aged_subset %>%
  filter(Prey_Com_Name == "Pacific Hake" & !is.na(Predator_ID)) %>%
  group_by(pred_ages, prey_ages) %>%
  summarize(prey_wt = sum(Prey_Weight_g)) %>%
  left_join(stomachs_n) %>%
  left_join(total_wt) %>%
  mutate(wt_prop = (prey_wt / total_wt))

# Calculate overall percentage by wt
mean(intrasp$wt_prop)

# Keep only the needed columns 
intrasp <- intrasp[, c("pred_ages", "prey_ages", "sample_size", "wt_prop")]

# Fill dataframe with missing predator and prey ages
all_ages <- data.frame(pred_ages = rep(1:15, each = 15), 
                       prey_ages = rep(1:15, 15), 
                       sample_size = rep(NA, 225), 
                       wt_prop = rep(NA, 225))

# Merge dataframe of all values w/ data, remove replicated rows, fill values
intrasp_full <- intrasp %>%
  full_join(all_ages) %>%
  arrange(pred_ages, prey_ages) %>%
  distinct(pred_ages, prey_ages, .keep_all = TRUE) %>%
  fill(sample_size)

# Replace remaining NAs with 0s
intrasp_full[is.na(intrasp_full)] <- 0

# Create new .csv with new values
# Update column names to match CEATTLE data input
colnames(intrasp_full)[c(3, 4)] <- c("Sample_size", "Stomach_proportion_by_weight")
write.csv(intrasp_full, "data/diet/full_hake_diet.csv", row.names = FALSE)

# Plot hake diet
df <- melt(intrasp[, -3], id.vars = c("pred_ages", "prey_ages"))

diet_plot <- ggplot(df, aes(x=as.factor(pred_ages), y=value, fill=as.factor(prey_ages))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits = factor(1:15)) +  # add in missing predator ages
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("predator hake age") + ylab("diet proportion by weight") +
  labs(fill = "prey hake age")
diet_plot

ggsave(filename = "plots/diet/cannibalism_overall.png", 
       diet_plot, width=200, height=120, units="mm", dpi=300)


### See if it's worth doing time-varying (yearly) predation -------------------
stomachs_yearly <- aged_subset %>%
  group_by(Year, pred_ages) %>%
  summarize(sample_size = n())

total_wt_yearly <- aged_subset %>%
  group_by(Year, pred_ages) %>%
  summarize(total_wt = sum(Prey_Weight_g, na.rm = TRUE))

intrasp_yearly <- aged_subset %>%
  filter(Prey_Com_Name == "Pacific Hake" & !is.na(Predator_ID)) %>%
  group_by(Year, pred_ages, prey_ages) %>%
  summarize(prey_wt = sum(Prey_Weight_g)) %>%
  left_join(stomachs_yearly) %>%
  left_join(total_wt_yearly) %>%
  mutate(wt_prop = (prey_wt / total_wt))

# Calculate overall percentage by wt to double check
mean(intrasp_yearly$wt_prop)

intrasp_yearly2 <- melt(intrasp_yearly[, c("Year", "pred_ages", "prey_ages", "wt_prop")],
                        id.vars = c("Year", "pred_ages", "prey_ages"))

diet_plot_yearly <- ggplot(intrasp_yearly2, aes(x=as.factor(pred_ages), y=value, fill=as.factor(prey_ages))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits = factor(1:15)) +  # add in missing predator ages
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("predator hake age") + ylab("diet proportion by weight") +
  labs(fill = "prey hake age") +
  facet_wrap(~Year)
diet_plot_yearly

ggsave(filename = "plots/diet/cannibalism_yearly.png", 
       diet_plot_yearly, width=300, height=200, units="mm", dpi=300)

