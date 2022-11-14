# Script for estimating hake age from length of predator and prey hake in the
# SWFSC dataset

library(ggplot2)
library(viridis)
library(FSA)
library(dplyr)
library(tidyr)
library(stats4)

# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

all_pred <- read.csv("data/diet/CCTD/hake_pred.csv")
all_prey <- read.csv("data/diet/CCTD/hake_prey.csv") 


### Parameterize length to age calculation  -----------------------------------
# Read in maturity data
maturity <- read.csv("~/Desktop/Local/hake-CEATTLE/Resources/hake-assessment-master/data/hake-maturity-data.csv")
age_length <- data.frame(na.omit(cbind(Age = maturity$Age, Length = maturity$Length_cm)))

hake_ages <- 0:15
hake_lengths <- min(age_length$Length):max(age_length$Length)

like <- function(logLinf,logK,a0,logSigma) {
  # Extract the parameters
  Linf <- exp(logLinf); K <- exp(logK); Sigma <- exp(logSigma)
  
  # make the model predictions
  Pred <- (-log(1 - Lengths/Linf) / K) + a0
  
  # Compute the negative log-likelihood
  NegLogL <- -1*sum(dnorm(Ages,Pred,Sigma,TRUE))
  
  return(NegLogL)
}  

start <- list(logLinf=log(100),logK=log(0.2),a0=0,logSigma=5)
fixed <- NULL
Ages <- age_length$Age
Lengths <- age_length$Length
mleOutput <- mle(like,start=start)
print(summary(mleOutput))

# Extract the parameters
Linf <- exp(coef(mleOutput)[1])
K <- exp(coef(mleOutput)[2])
a0 <- coef(mleOutput)[3]
Sigma <- exp(coef(mleOutput)[4])

# Check AIC... yike
AIC(mleOutput)

# Plot the model fit
# Generate modelled data
predicted_ages <- (-log(1 - hake_lengths/Linf) / K) + a0
predicted_data <- data.frame(hake_lengths, predicted_ages)

ggplot(predicted_data, aes(x = hake_lengths, y = predicted_ages)) +
  geom_line() +
  geom_point(data = age_length, aes(x = Lengths, y = Age))


### Calculate predator ages ---------------------------------------------------
pred_ages <- (-log(1 - all_pred$FL_cm / Linf) / K) + a0

max(na.omit(pred_ages))  # check maximum
min(na.omit(pred_ages))  # check minimum
pred_ages[pred_ages < 1] <- 1  # replace values < 1 (lower accumulation age) <---- THIS SHOULD CHANGE WHEN AGE 0 WORKS
pred_ages[pred_ages > 15] <- 15  # set upper accumulation age to 15
pred_ages <- round(pred_ages, digits = 0)  # round to whole number

# Add ages column to predator dataset 
new_pred <- cbind(all_pred, pred_ages)


### Run prey age calculation --------------------------------------------------
prey_ages <- (-log(1 - (all_prey$Prey_Length1 / 10) / Linf) / K) + a0

max(na.omit(prey_ages))  # check maximum
min(na.omit(prey_ages))  # check minimum
prey_ages[prey_ages < 1] <- 1  # replace values < 1 (lower accumulation age) <---- THIS SHOULD CHANGE WHEN AGE 0 WORKS
prey_ages <- round(prey_ages, digits = 0)  # round to whole number

# Add ages column to prey dataset - includes non-hake prey
new_prey <- cbind(all_prey, prey_ages)
# Filter prey dataset to only include hake
new_hake_prey <- new_prey %>% filter(Prey_Com_Name == "Pacific Hake")


### Read in FEAT data and compare fit -----------------------------------------
FEAT_data <- read.csv("data/diet/FEAT/all_FEAT_diet.csv")
FEAT_prey_lengths <- read.csv("data/diet/FEAT/hake_prey_lengths_FEAT.csv")

# Add prey lengths based on stomach ID
FEAT_hake <- merge(FEAT_data, FEAT_prey_lengths, by = "stomach_uuid", all.x = TRUE) %>%
  filter(prey_category_long == "Gadiformes") %>%
  select(predator_length_cm, predator_age, measure_value)

FEAT_hake$prey_ages <- (-log(1 - FEAT_hake$measure_value/Linf) / K) + a0
FEAT_hake$prey_ages[FEAT_hake$prey_ages < 1] <- 1  # replace values < 1 (lower accumulation age) <---- THIS SHOULD CHANGE WHEN AGE 0 WORKS

FEAT_ages <- as.data.frame(rbind(cbind(age = FEAT_hake$predator_age, length = FEAT_hake$predator_length_cm,
                                       data = rep("FEAT predators", length(FEAT_hake$predator_ages)))))

# Plot fit --------------------------------------------------------------------
all_ages <- as.data.frame(rbind(cbind(age = maturity$Age,
                                      length = maturity$Length_cm,
                                      data = rep("assessment", length(maturity$Age))),
                                cbind(age = new_pred$pred_ages,
                                      length = new_pred$FL_cm,
                                      data = rep("CCTD predators", length(new_pred$pred_ages))),
                                cbind(age = new_hake_prey$prey_ages,
                                      length = (new_hake_prey$Prey_Length1 / 10),  # prey are in mm
                                      data = rep("CCTD prey", length(new_hake_prey$prey_ages))),
                                cbind(age = FEAT_hake$predator_age,
                                      length = FEAT_hake$predator_length_cm,
                                      data = rep("FEAT predators", length(FEAT_hake$predator_age))),
                                cbind(age = FEAT_hake$prey_ages,
                                      length = FEAT_hake$measure_value,
                                      data = rep("FEAT prey", length(FEAT_hake$prey_ages)))))
all_ages$age <- as.numeric(all_ages$age)
all_ages$length <- as.numeric(all_ages$length)

all_ages <- na.omit(all_ages)

growth_curve <- ggplot(all_ages, aes(x = age, y = length, color = data, shape = data)) +
  geom_point(alpha = 0.4, size = 3) +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9)
growth_curve

ggsave(filename = "plots/diet/growth_curve.png", growth_curve,
       bg = "transparent", width=160, height=80, units="mm", dpi=300)


# Write aged datasets to file -------------------------------------------------
write.csv(new_pred, "data/diet/CCTD/hake_aged_pred.csv", row.names = FALSE)
write.csv(new_prey, "data/diet/CCTD/hake_aged_prey.csv", row.names = FALSE)


# ### Age-length key method for ageing ------------------------------------------
# # Set up hake age-length key
# hake_lbin <- c(0, seq(20, 68, by = 2), 999)
# hake_age_bin <- c(0:14 + 0.5, 99)
# hake_ages <- 1:15
# 
# maturity$BIN <- cut(maturity$Length_cm, breaks = hake_lbin)
# levels(maturity$BIN) <- 1:length(hake_lbin[-1])
# 
# maturity$AgeBIN <- cut(maturity$Age, breaks = hake_age_bin)
# levels(maturity$AgeBIN) <- hake_ages
# 
# hake <- maturity[-which(is.na(maturity$AgeBIN)),]
# hake_table <- with(hake,table(BIN,AgeBIN))
# alk_hake <- prop.table(hake_table, margin=1)
# alk_hake[is.na(alk_hake)] <- 0
# 
# # # Save for age_trans_matrix in CEATTLE input excel sheet
# age_trans_matrix <- t(as.data.frame.matrix(alk_hake))
# write.csv(age_trans_matrix, "data/assessment/age_trans_matrix.csv", row.names = FALSE)
# 
# # Predator age calculations
# pred_ages <- alkIndivAge(key = alk_hake,
#                          formula = age ~ length,
#                          data = data.frame(cbind(length = all_pred$FL_cm, age = rep(NA))),
#                          type="CR")
# # Check min and max age
# max(pred_ages$age)
# min(pred_ages)
# 
# # Add to exising data
# all_pred$age <- pred_ages$age
# 
# # Prey age calculations - more complicated
# prey_hake <- all_prey %>%
#   filter(Prey_Com_Name == "Pacific Hake")
# prey_hake$Prey_Length1 <- prey_hake$Prey_Length1 / 10  # convert to cm
# prey_hake[is.na(prey_hake)] <- 1  # All NAs in length correspond to immature - below age 2
# 
# prey_ages <- alkIndivAge(key = alk_hake,
#                          formula = age ~ length,
#                          data = data.frame(cbind(length = prey_hake$Prey_Length1, age = rep(NA))),
#                          type="CR")
# # Check min and max age
# max(prey_ages$age)
# min(prey_ages$age)
# 
# # Add to existing data
# prey_hake$age <- prey_ages$age
# all_prey$age <- NA
# all_prey$age[all_prey$Prey_Com_Name == "Pacific Hake"] <- prey_ages$age