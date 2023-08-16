# Script for estimating hake age from length of predator and prey hake in the
# SWFSC dataset

library(ggplot2)
library(viridis)
library(FSA)
library(dplyr)
library(tidyr)
library(stats4)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

all_pred <- read.csv("data/diet/CCTD/hake_pred.csv")
all_prey <- read.csv("data/diet/CCTD/hake_prey.csv") 

# Read in age data from the survey & filter to only be hake
survey_ages <- read.csv("data/diet/acoustic_survey.csv") %>% 
  filter(common_name == "Pacific hake")

age_length <- na.omit(cbind.data.frame(survey = survey_ages$survey,
                                       age = survey_ages$age, 
                                       length = survey_ages$length)) %>%
  mutate(Year = substr(survey, 1, 4))

# Plot age/length data by year
growth_yearly <- ggplot(age_length, aes(x = age, y = length)) +
  geom_point() +
  ylab("Length (cm)") +
  facet_wrap(~Year, ncol = 6)
growth_yearly
# ggsave("plots/diet/growth_yearly.png", growth_yearly, 
#        width=250, height = 120, units = "mm", dpi=300)


### Try different growth curve calculations -----------------------------------
# Plot von Bertalanffy growth curve 
hake_ages <- 0:20
hake_lengths <- min(age_length$length):max(age_length$length)

like <- function(logLinf, logK, a0, logSigma) {
  # Extract the parameters
  Linf <- exp(logLinf); K <- exp(logK); Sigma <- exp(logSigma)
  
  # make the model predictions
  Pred <- (-log(1 - Lengths/Linf) / K) + a0
  
  # Compute the negative log-likelihood
  NegLogL <- -1*sum(dnorm(Ages, Pred, Sigma, TRUE))
  
  return(NegLogL)
}  

# FSA::vbStarts(length ~ age, data = age_length)  # suggested starting params
start <- list(logLinf = log(90), logK=log(0.4), a0=0, logSigma=5)
fixed <- NULL
Ages <- age_length$age
Lengths <- age_length$length
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

ggplot(predicted_data, aes(x = predicted_ages, y = hake_lengths)) +
  geom_line() +
  geom_point(data = age_length, aes(x = age, y = length))


# Try Schnute parameterization w/ both normal & lognormal VBGF ----------------
a1 <- 1
a2 <- 20

FSA::vbStarts(length ~ age, data = age_length, param = "Schnute", ages2use = c(a1, a2))

# Normal
vbgf.nls2 <- nls(length ~ la1 + (la2 - la1) * 
                   (1 - exp(-k * (age - a1))) / (1 - exp(-k * (a2 - a1))), 
                 data = age_length, 
                 start = list(la1 = 23, 
                              la2 = 56, 
                              k = 0.4))
summary(vbgf.nls2)
summary(vbgf.nls2, cor=TRUE)$correlation %>% round(2)

plot(fitted(vbgf.nls2), resid(vbgf.nls2))
abline(h=0)

nls <- coef(vbgf.nls2)

ggplot(age_length) +
  geom_point(aes(x = age, y = length), alpha = 0.25) +
  geom_line(aes(x = age, y = nls[1] + (nls[2] - nls[1]) * (1-exp(-nls[3] * (age - 1))) /
                  (1-exp(-nls[3]*14))))

k <- nls[3]
la1 <- nls[1]
la2 <- nls[2]

# Lognormal
vbgf.loglike <- function(log.pars, dat, a1, a2) {
  pars <- exp(log.pars)
  
  l.pred <- pars["la1"] + (pars["la2"] - pars["la1"]) * 
    (1-exp(-pars["k"]*(dat$age - a1))) / 
    (1-exp(-pars["k"]*(a2-a1)))
  
  nll <- -dlnorm(x = dat$length, 
                 meanlog = log(l.pred) - pars['cv']^2, 
                 sdlog = pars['cv'], 
                 log = TRUE) %>%
    sum()
  return(nll)
}

dat <- filter(age_length, !is.na(age))
pars.init <- log(c(la1 = 24, la2 = 56, k = 0.36, cv = 0.4))
vbgf.optim <- optim(pars.init, vbgf.loglike, dat=dat, a1=a1, a2=a2)

exp(vbgf.optim$par)[1:4]

nls <- coef(vbgf.nls2)
optim <- exp(vbgf.optim$par)

# Solve for Linf and a0 using Schnute parameters
Linf2 <- la2 - la1 * exp(-k * (a2 - a1)) / (1 - exp(-k * (a2 - a1)))
a0_2 <- a1 + (1/k * log((la2 - la1) / (la2 - (la1 * exp(-k * (a2 - a1))))))
 
predicted_ages2 <- (-log(1 - hake_lengths / Linf2) / k) + a0_2
predicted_data2 <- data.frame(hake_lengths, predicted_ages2)

ggplot() +
  geom_point(data = age_length, aes(x = age, y = length), color = "gray", alpha = 0.5) +
  geom_line(data = age_length, aes(x = age, y = nls[1] + (nls[2] - nls[1]) * 
                              (1-exp(-nls[3]*(age-1))) / 
                              (1-exp(-nls[3]*14)), 
                            linetype = "Schnute - normal", color = "Schnute - normal"),
            linewidth = 1.5) +
  geom_line(data = age_length, aes(x = age, y = optim[1] + (optim[2] - optim[1]) *
                              (1-exp(-optim[3]*(age-1))) / 
                              (1-exp(-optim[3]*14)), 
                            linetype = "Schnute - lognormal", color = "Schnute - lognormal"),
            linewidth = 1.5) +
  geom_line(data = predicted_data, aes(x = predicted_ages, y = hake_lengths,
                                linetype = "VBGF", color = "VBGF"),
            linewidth = 1.5) +
  geom_line(data = predicted_data2, aes(x = predicted_ages2, y = hake_lengths,
                                       linetype = "VBGF w/ Schnute", color = "VBGF w/ Schnute"),
            linewidth = 1.5) +
  ylab("Length (cm)") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  labs(linetype = "Model")


# # Find von Bert parameters for specific years -------------------------------
# VBGF_params <- function(year) {
#   df <- age_length %>% filter(Year == year)
#   hake_ages <- 0:20
#   hake_lengths <- min(df$length):max(df$length)
#   
#   like <- function(logLinf, logK, a0, logSigma) {
#     # Extract the parameters
#     Linf <- exp(logLinf); K <- exp(logK); Sigma <- exp(logSigma)
#     
#     # make the model predictions
#     Pred <- (-log(1 - Lengths/Linf) / K) + a0
#     
#     # Compute the negative log-likelihood
#     NegLogL <- -1*sum(dnorm(Ages, Pred, Sigma, TRUE))
#     
#     return(NegLogL)
#   }  
#   
#   start_est <- FSA::vbStarts(length ~ age, data = df)  # suggested starting params
#   start <- list(logLinf = log(90), logK=log(0.4), a0=0, logSigma=6)
#   Ages <- df$age
#   Lengths <- df$length
#   mleOutput <- mle(like, start=start)
#   
#   # Extract the parameters
#   Linf <- exp(coef(mleOutput)[1])
#   K <- exp(coef(mleOutput)[2])
#   a0 <- coef(mleOutput)[3]
#   Sigma <- exp(coef(mleOutput)[4])
#   
#   # Create dataframe of parameters
#   parameters <- cbind.data.frame(
#     year = year,
#     Linf = exp(coef(mleOutput)[1]),
#     K = exp(coef(mleOutput)[2]),
#     a0 = coef(mleOutput)[3],
#     Sigma = Sigma <- exp(coef(mleOutput)[4])
#   )
#   
#   # Plot the model fit
#   predicted_ages <- (-log(1 - hake_lengths/Linf) / K) + a0
#   predicted_data <- data.frame(hake_lengths, predicted_ages)
#   
#   plot <- ggplot(predicted_data, aes(x = predicted_ages, y = hake_lengths)) +
#     geom_line() +
#     geom_point(data = age_length, aes(x = age, y = length)) +
#     ggtitle(year)
#   
#   return(list(params = parameters, plot = plot))
# }
# 
# years_used <- c(1989, 1992, 1995, 1998, 2005, 2007, 2011, 2015, 2017)
# test <- VBGF_params(1992)
# test$plot


### Calculate predator ages ---------------------------------------------------
pred_ages <- (-log(1 - all_pred$FL_cm / Linf) / K) + a0

max(na.omit(pred_ages))  # check maximum
min(na.omit(pred_ages))  # check minimum
pred_ages[pred_ages < 1] <- 1  # replace values < 1 (lower accumulation age)
pred_ages[pred_ages > 20] <- 20  # set upper accumulation age to 20
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
all_ages <- as.data.frame(rbind(cbind(age = age_length$age,
                                      length = age_length$length,
                                      data = rep("assessment", length(age_length$age))),
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

growth_curve <- ggplot(all_ages, aes(x = age, y = length, 
                                     color = data, shape = data, alpha = data)) +
  geom_point(size = 3) +
  scale_alpha_discrete(range = c(0.2, 1)) +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9)
growth_curve

ggsave(filename = "plots/diet/growth_curve.png", growth_curve,
       bg = "transparent", width=160, height=80, units="mm", dpi=300)


# Write aged datasets to file -------------------------------------------------
write.csv(new_pred, "data/diet/CCTD/hake_aged_pred.csv", row.names = FALSE)
write.csv(new_prey, "data/diet/CCTD/hake_aged_prey.csv", row.names = FALSE)


### Age-length key method for ageing ------------------------------------------
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
# # Save for age_trans_matrix in CEATTLE input excel sheet
# age_trans_matrix <- t(as.data.frame.matrix(alk_hake))
# # write.csv(age_trans_matrix, "data/assessment/age_trans_matrix.csv", row.names = FALSE)
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
# # Add to existing data
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
