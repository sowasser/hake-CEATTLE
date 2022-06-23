# Script for estimating hake age from length of predator and prey hake in the
# SWFSC dataset

library(ggplot2)
library(ggsidekick)
library(viridis)
library(FSA)
library(dplyr)
library(tidyr)

all_pred <- read.csv("data/diet/Full dataset/hake_pred.csv")
all_prey <- read.csv("data/diet/Full dataset/hake_prey.csv")

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

# Add ages column to prey dataset - includes non-hake prey
new_prey <- cbind(all_prey, prey_ages)
# Round to whole number 
new_prey$prey_ages <- round(new_prey$prey_ages, digits = 0)
# Filter prey dataset to only include hake
new_hake_prey <- new_prey %>% filter(Prey_Com_Name == "Pacific Hake")

# Plot fit --------------------------------------------------------------------
all_ages <- as.data.frame(rbind(cbind(age = maturity$Age, length = maturity$Length_cm, 
                                      data = rep("original", length(maturity$Age))),
                                cbind(age = new_pred$pred_ages, length = new_pred$FL_cm, 
                                      data = rep("predator hake", length(new_pred$pred_ages))),
                                cbind(age = new_hake_prey$prey_ages, length = (new_hake_prey$Prey_Length1 / 10),  # prey are in mm
                                      data = rep("prey hake", length(new_hake_prey$prey_ages)))))
all_ages$age <- as.numeric(all_ages$age)
all_ages$length <- as.numeric(all_ages$length)

all_ages <- na.omit(all_ages)

growth_curve <- ggplot(all_ages, aes(x = age, y = length, color = data)) +
  geom_point(alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek()
growth_curve

ggsave(filename = "plots/diet/growth_curve.png", growth_curve, 
       width=200, height=120, units="mm", dpi=300)


# Write aged datasets to file -------------------------------------------------
write.csv(new_pred, "data/diet/Full dataset/hake_aged_pred.csv", row.names = FALSE)
write.csv(new_prey, "data/diet/Full dataset/hake_aged_prey.csv", row.names = FALSE)
