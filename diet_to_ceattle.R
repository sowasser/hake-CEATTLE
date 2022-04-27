# Script for updating the hake intraspecies predation data for use in CEATTLE

library(r4ss)
library(readxl)
library(dplyr)
library(tidyr)

all_hake <- read.csv("data/diet/full_hake_diet.csv")
all_pred <- read.csv("data/diet/full_hake_pred.csv")
all_prey <- read.csv("data/diet/full_prey.csv")

### Create weight-age keys ----------------------------------------------------
hake_ages <- 0:15

# Mean weights from empirical weight-at-age table
hake_wt <- c(0.022340741,	0.111942593,	0.258737037,	0.381925926,	0.472466667,	
             0.531718519,	0.58675,	0.647811111,	0.710514815,	0.765505556,	
             0.824953704,	0.928575926,	0.986738889,	1.024092593,
             1.082364815,	1.212831481, 10)

# Plot just to make sure everything makes sense
test <- as.data.frame(cbind(hake_ages, hake_wt))
ggplot(test, aes(x = hake_ages, y = hake_wt)) +
  geom_point() +
  theme_sleek()

# Cut weight data into age categories
pred_ages <- cut(all_pred$Predator_Weight_kg, breaks = hake_wt, right = TRUE, labels = hake_ages)

new_pred <- cbind(all_pred, pred_ages)
