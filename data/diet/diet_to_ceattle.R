# Script for updating the hake intraspecies predation data for use in CEATTLE

library(ggplot2)
library(ggsidekick)
library(FSA)
library(dplyr)
library(tidyr)

all_pred <- read.csv("data/diet/Full dataset/full_hake_pred.csv")
all_prey <- read.csv("data/diet/Full dataset/full_prey.csv")

### Parameterize length to age calculation  -----------------------------------
hake_ages <- 0:15

# Read in maturity data
maturity <- read.csv("~/Desktop/Local/hake-CEATTLE/Resources/hake-assessment-master/data/hake-maturity-data.csv")

# Estimate VBGF (http://derekogle.com/fishR/2019-12-31-ggplot-vonB-fitPlot-1)
vb <- vbFuns(param="Typical")  # define von Bert function
# Get reasonable starting values, fit model, extract parameters
f.starts <- vbStarts(Length_cm ~ Age, data = maturity)
f.fit <- nls(Length_cm ~ vb(Age, Linf, K, t0), data = maturity, start = f.starts) 
params <- coef(f.fit) 

# Calculate ages from lengths in dataset
age_calc <- function(lengths, Linf, K, t0) {
  ages <- c()
  for(L in lengths) {
    a <- max(1, ((-log(1 - L/Linf))/K + t0))
    ages <- c(ages, a)
  }
  
  return(ages)
}

### Run predator age calculation ----------------------------------------------
pred_ages <- age_calc(lengths = all_pred$FL_cm, 
                      Linf = params[1], K = params[2], t0 = params[3])

# Add ages column to predator dataset 
new_pred <- cbind(all_pred, pred_ages)
# Fill in age = 15 for any length > Linf
new_pred$pred_ages[new_pred$FL_cm > params[1]] <- 15
# Round to whole number 
new_pred$pred_ages <- round(new_pred$pred_ages, digits = 0)
# Set any ages > 15 to 15 (accumulator age)
new_pred$pred_ages[new_pred$pred_ages > 15] <- 15

### Run prey age calculation --------------------------------------------------
prey_ages <- age_calc(lengths = (all_prey$Prey_Length1 / 10),  # prey are in mm
                      Linf = params[1], K = params[2], t0 = params[3])

# Add ages column to prey dataset 
new_prey <- cbind(all_prey, prey_ages)
# Round to whole number 
new_prey$prey_ages <- round(new_prey$prey_ages, digits = 0)
# Replace any non-hake prey items with NA
new_prey$prey_ages[new_prey$Prey_Com_Name != "Pacific Hake"] <- NA


### Combine into new dataset & prep for use in CEATTLE ------------------------
aged_dataset <- merge(new_pred, new_prey, all = TRUE)

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

# Get number of stomachs per predator age
stomachs_n <- aged_subset %>%
  group_by(pred_ages) %>%
  summarize(sample_size = n())

cannibalism <- aged_subset %>%
  filter(Prey_Com_Name == "Pacific Hake")

# Get proportion by weight for each predator & prey age
diet_prop_overall <- aged_subset %>% 
  group_by(pred_ages, prey_ages) %>%
  summarize(sample_size = n()) %>%
  mutate(prop = sample_size / sum(sample_size))

