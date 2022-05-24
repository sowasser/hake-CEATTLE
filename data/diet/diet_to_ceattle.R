# Script for updating the hake intraspecies predation data for use in CEATTLE

library(ggplot2)
library(ggsidekick)
library(viridis)
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


# Plot fit --------------------------------------------------------------------
all_ages <- as.data.frame(rbind(cbind(age = maturity$Age, length = maturity$Length_cm, 
                                      data = rep("original", length(maturity$Age))),
                                cbind(age = new_pred$pred_ages, length = new_pred$FL_cm, 
                                      data = rep("predator hake", length(new_pred$pred_ages))),
                                cbind(age = new_prey$prey_ages, length = new_prey$Prey_Length1, 
                                      data = rep("prey_hake", length(new_prey$prey_ages)))))
all_ages$age <- as.numeric(all_ages$age)
all_ages$length <- as.numeric(all_ages$length)

all_ages <- na.omit(all_ages)

growth_curve <- ggplot(all_ages, aes(x = age, y = length, color = data)) +
  geom_point(alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek()
growth_curve



### Combine into new dataset --------------------------------------------------
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

