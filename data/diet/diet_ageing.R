# Script for estimating hake age from length of predator and prey hake in the
# SWFSC dataset

library(ggplot2)
library(viridis)
library(FSA)
library(dplyr)
library(tidyr)

# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

all_pred <- read.csv("data/diet/CCTD/hake_pred.csv")
all_prey <- read.csv("data/diet/CCTD/hake_prey.csv") 


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


# Read in FEAT data -----------------------------------------------------------
FEAT_data <- read.csv("data/diet/FEAT/all_FEAT_diet.csv")
FEAT_prey_lengths <- read.csv("data/diet/FEAT/hake_prey_lengths_FEAT.csv")

# Add prey lengths based on stomach ID
FEAT_hake <- merge(FEAT_data, FEAT_prey_lengths, by = "stomach_uuid", all.x = TRUE) %>%
  filter(prey_category_long == "Gadiformes") %>%
  select(predator_length_cm, predator_age, measure_value)

FEAT_hake$prey_ages <- age_calc(lengths = FEAT_hake$measure_value,
                                Linf = params[1], K = params[2], t0 = params[3])

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
  geom_point(alpha = 0.3, size = 3) +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9)
growth_curve

ggsave(filename = "plots/diet/growth_curve.png", growth_curve,
       bg = "transparent", width=160, height=80, units="mm", dpi=300)


# Write aged datasets to file -------------------------------------------------
write.csv(new_pred, "data/diet/CCTD/hake_aged_pred.csv", row.names = FALSE)
write.csv(new_prey, "data/diet/CCTD/hake_aged_prey.csv", row.names = FALSE)


# ### Age-length key method for ageing ------------------------------------------
# # Mean Length at age (this is Mn_LatAge in the CEATTLE input excel sheet)
# mean_lengths <- maturity %>% group_by(Age) %>% 
#   summarise(Length = mean(Length_cm))
# 
# # Set up hake age-length key
# hake_lbin <- c(0, seq(23, 56, by = 1), 999)
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
# # Predator age calculations
# pred_ages <- alkIndivAge(key = alk_hake, 
#                          formula = age ~ length, 
#                          data = data.frame(cbind(length = all_pred$FL_cm, age = rep(NA))), 
#                          type="CR")
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
