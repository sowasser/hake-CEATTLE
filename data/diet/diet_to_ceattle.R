# Script for updating the hake intraspecies predation data for use in CEATTLE

library(tidyverse)
library(ggsidekick)
library(viridis)
library(reshape2)

# Read in full aged dataset
aged_dataset <- read.csv("data/diet/CCTD_FEAT_combined.csv")

# Look at instances where prey hake age = NA - all are immature
CCTD_prey <- read.csv("data/diet/CCTD/hake_aged_prey.csv")
  
CCTD_prey %>%
  filter(Prey_Com_Name == "Pacific Hake") %>%
  filter(is.na(prey_ages))

# Replace those NAs with age 1
aged_dataset$prey_age[is.na(aged_dataset$prey_age) & aged_dataset$prey_name == "Pacific Hake"] <- 1

# Get total stomach weights for each predator and proportional weight for prey hake
aged_wt <- aged_dataset[, c("Predator_ID", "year", "predator_age", "prey_name", "prey_age", "prey_wt")] %>%
  group_by(Predator_ID) %>%
  mutate(stomach_wt = sum(prey_wt, na.rm = TRUE)) %>%
  filter(prey_name == "Pacific Hake") %>%
  mutate(hake_prey_prop = prey_wt / stomach_wt) %>%
  ungroup() %>%
  distinct()  # Remove duplicate rows - same pred ID, multiple hake prey

aged_wt$new_ID <- c(1:nrow(aged_wt))

# Remove duplicate rows - same pred ID, multiple hake prey but same prey wts - generalized for that prey species
hake_prop_wide <- aged_wt %>%
  select(new_ID, predator_age, prey_age, hake_prey_prop) %>%
  distinct() %>%
  pivot_wider(id_cols = c(new_ID, predator_age),  # columns to stay the same
              names_from = prey_age,  # columns to convert to wide
              values_from = hake_prey_prop,  # values to fill in
              values_fill = 0)  # what to fill for missing values

# Re-order and re-name columns
hake_prop_wide <- hake_prop_wide[, c(2, 4, 5, 6, 3)]
colnames(hake_prop_wide)[2:5] <- c("prey_a1", "prey_a2", "prey_a3", "prey_a5")


# Create overall intraspecies predation dataset -------------------------------
# Get number of stomachs per predator age
stomachs_n <- aged_dataset %>%
  group_by(predator_age) %>%
  summarize(sample_size = n())
  
# Total stomach weight per predator age
total_wt <- aged_dataset %>%
  group_by(predator_age) %>%
  summarize(total_wt = sum(prey_wt, na.rm = TRUE))

# Combine summarized datasets and calculate stomach weight proportions (overall)
intrasp <- aged_dataset %>%
  filter(prey_name == "Pacific Hake" & !is.na(Predator_ID)) %>%
  group_by(predator_age, prey_age) %>%
  summarize(prey_wt = sum(prey_wt)) %>%
  left_join(stomachs_n) %>%
  left_join(total_wt) %>%
  mutate(wt_prop = (prey_wt / total_wt))

# Calculate overall percentage by wt
mean(intrasp$wt_prop)

# Keep only the needed columns 
intrasp <- intrasp[, c("predator_age", "prey_age", "sample_size", "wt_prop")]

# Fill dataframe with missing predator and prey ages
all_ages <- data.frame(predator_age = rep(1:15, each = 15), 
                       prey_age = rep(1:15, 15), 
                       sample_size = rep(NA, 225), 
                       wt_prop = rep(NA, 225))

# Merge dataframe of all values w/ data, remove replicated rows, fill values
intrasp$predator_age <- as.integer(intrasp$predator_age)
intrasp_full <- intrasp %>%
  full_join(all_ages) %>%
  arrange(predator_age, prey_age) %>%
  distinct(predator_age, prey_age, .keep_all = TRUE) %>% 
  fill(sample_size)  # had an issue with this line, restarted R & it works!

# Replace remaining NAs with 0s
intrasp_full[is.na(intrasp_full)] <- 0

# Plot hake diet
df <- melt(intrasp[, -3], id.vars = c("predator_age", "prey_age"))

diet_plot <- ggplot(df, aes(x=as.factor(predator_age), y=value, fill=as.factor(prey_age))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits = factor(1:15)) +  # add in missing predator ages
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("predator hake age") + ylab("diet proportion by weight") +
  labs(fill = "prey hake age")
diet_plot

ggsave(filename = "plots/diet/cannibalism_overall.png", 
       diet_plot, width=200, height=120, units="mm", dpi=300)


### Reshape data for Dirichlet modeling ---------------------------------------
for_dirichlet <- aged_dataset %>%
  select(Predator_ID, year, predator_age, prey_name, prey_age, prey_wt) %>%
  # Find total stomach weight for each predator
  group_by(Predator_ID) %>%
  mutate(stomach_wt = sum(prey_wt, na.rm = TRUE)) %>% 
  # Select only prey hake and find weight proportion
  filter(prey_name == "Pacific Hake") %>%  
  mutate(hake_prey_prop = prey_wt / stomach_wt) %>%
  ungroup() %>%
  distinct()  %>%  # Remove duplicates - keep only one record per predator
  mutate(new_ID = c(1:length(Predator_ID))) %>%
  select(new_ID, predator_age, prey_age, hake_prey_prop) %>%
  pivot_wider(id_cols = c(new_ID, predator_age),  # rows to stay the same
              names_from = prey_age,  # columns to convert to wide
              values_from = hake_prey_prop,  # values to fill in
              values_fill = 0) %>%  # what to fill for missing values 
  mutate(predator_age = as.integer(predator_age)) %>%
  arrange(predator_age)

# Re-order and rename data for final dataset
for_dirichlet <- for_dirichlet[, c(2, 4, 5, 6, 3)]
colnames(for_dirichlet)[2:5] <- c("prey_a1", "prey_a2", "prey_a3", "prey_a5")
unique(as.character(for_dirichlet$predator_age))
for_dirichlet$predator_age <- recode(for_dirichlet$predator_age, 
                                   "1" = "pred_a1", "2" = "pred_a2",
                                   "3" = "pred_a3", "4" = "pred_a4",
                                   "5" = "pred_a5", "6" = "pred_a6",
                                   "7" = "pred_a7", "8" = "pred_a8",
                                   "9" = "pred_a9", "10" = "pred_a10",
                                   "11" = "pred_a11", "12" = "pred_a12",
                                   "13" = "pred_a13", "14" = "pred_a14",
                                   "15" = "pred_a15")

# Add "wtg" column from Dirichlet example dataset
for_dirichlet <- cbind(for_dirichlet[, 1], 
                        wtg = rep(1, length(for_dirichlet[, 1])),
                        for_dirichlet[, 2:5])

# Add column with remaining diet proportion - assigned to "other" prey
for_dirichlet$other <- 1 - rowSums(for_dirichlet[, 3:6])

colnames(for_dirichlet)[1] <- "Predator"
write.csv(for_dirichlet, "data/diet/Dirichlet/hake_for_dirichlet.csv", row.names = FALSE)
  

### See if it's worth doing time-varying (yearly) predation -------------------
stomachs_yearly <- aged_dataset %>%
  group_by(year, predator_age) %>%
  summarize(sample_size = n())

total_wt_yearly <- aged_dataset %>%
  group_by(year, predator_age) %>%
  summarize(total_wt = sum(prey_wt, na.rm = TRUE))

intrasp_yearly <- aged_dataset %>%
  filter(prey_name == "Pacific Hake" & !is.na(Predator_ID)) %>%
  group_by(year, predator_age, prey_age) %>%
  summarize(prey_wt = sum(prey_wt)) %>%
  left_join(stomachs_yearly) %>%
  left_join(total_wt_yearly) %>%
  mutate(wt_prop = (prey_wt / total_wt))

# Calculate overall percentage by wt to double check
mean(intrasp_yearly$wt_prop)

intrasp_yearly <- melt(intrasp_yearly[, c("year", "predator_age", "prey_age", "wt_prop")],
                       id.vars = c("year", "predator_age", "prey_age"))

diet_plot_yearly <- ggplot(intrasp_yearly, aes(x=as.factor(predator_age), y=value, fill=as.factor(prey_age))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits = factor(1:15)) +  # add in missing predator ages
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("predator hake age") + ylab("diet proportion by weight") +
  labs(fill = "prey hake age") +
  facet_wrap(~year)
diet_plot_yearly

ggsave(filename = "plots/diet/cannibalism_yearly.png", 
       diet_plot_yearly, width=300, height=200, units="mm", dpi=300)


### Get data ready to be added directly to CEATTLE ----------------------------
intrasp_ceattle <- cbind(Pred = rep(1, nrow(intrasp_full)),
                         Prey = rep(1, nrow(intrasp_full)),
                         Pred_sex = rep(0, nrow(intrasp_full)),
                         Prey_sex = rep(0, nrow(intrasp_full)),
                         Pred_age = intrasp_full$predator_age,
                         Prey_age = intrasp_full$prey_age,
                         Year = rep(0, nrow(intrasp_full)),
                         Sample_size = intrasp_full$sample_size,
                         Stomach_proportion_by_weight = intrasp_full$wt_prop)

write.csv(intrasp_ceattle, "data/diet/diet_for_CEATTLE.csv", row.names = FALSE)


### Read in & update post-Dirichlet dataset -----------------------------------
post_dirichlet <- read.csv("data/diet/Dirichlet/Dirichlet_Results.csv")

# Fix issues from export from the Dirichlet script
colnames(post_dirichlet) <- c("predator", "prey", "lower95", "upper95", "mode", 
                              "simple_average", "boot_average", "boot_mode",
                              "normalized_mode")
post_dirichlet <- post_dirichlet[-1, ]

# Remove "other" data column from analysis
post_dirichlet <- post_dirichlet %>% filter(prey != "other")
post_dirichlet[, 3:9] <- sapply(post_dirichlet[, 3:9], as.numeric)  # make numeric

# Plot different statistics from the analysis
post_dirichlet_long <- melt(post_dirichlet, id.vars = c("predator", "prey", "lower95", "upper95"))

dirichlet_results <- ggplot(post_dirichlet_long, aes(x = prey)) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), 
                width = .3, position = position_dodge(.9)) +
  geom_point(aes(x = prey, y = value, color = variable, shape = variable), size = 7, alpha = 0.5) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab("") + ylab("stomach proportion") +
  facet_wrap(~ predator)
dirichlet_results

ggsave(filename = "plots/diet/Dirichlet_results.png", 
       dirichlet_results, width=300, height=200, units="mm", dpi=300)

# Re-organize dataframe for CEATTLE
dirichlet_ceattle <- data.frame(cbind(Pred = rep(1, nrow(post_dirichlet)),
                                      Prey = rep(1, nrow(post_dirichlet)),
                                      Pred_sex = rep(0, nrow(post_dirichlet)),
                                      Prey_sex = rep(0, nrow(post_dirichlet)),
                                      Pred_age = post_dirichlet$predator,
                                      Prey_age = post_dirichlet$prey,
                                      Year = rep(0, nrow(post_dirichlet)),
                                      Sample_size = rep(10, nrow(post_dirichlet)),
                                      Stomach_proportion_by_weight = post_dirichlet$boot_average))

write.csv(dirichlet_ceattle, "data/diet/Dirichlet_for_CEATTLE.csv", row.names = FALSE)

