# Script for updating the hake intraspecies predation data for use in CEATTLE

library(ggplot2)
library(ggsidekick)
library(FSA)
library(dplyr)

all_hake <- read.csv("data/diet/full_hake_diet.csv")
all_pred <- read.csv("data/diet/full_hake_pred.csv")
all_prey <- read.csv("data/diet/full_prey.csv")

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

### Run prey age calculation --------------------------------------------------
prey_ages <- age_calc(lengths = (all_prey$Prey_Length1 / 10),  # prey are in mm
                      Linf = params[1], K = params[2], t0 = params[3])

# Add ages column to predator dataset 
new_prey <- cbind(all_prey, prey_ages)
# Fill in age = 15 for any length > Linf
new_prey$prey_ages[new_prey$FL_cm > params[1]] <- 15
# Round to whole number 
new_prey$prey_ages <- round(new_prey$prey_ages, digits = 0)

### Combine into new dataset & prep for use in CEATTLE ------------------------
aged_datset <- all_hake <- merge(new_pred, new_prey, all = TRUE)
