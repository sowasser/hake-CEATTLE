# Which predators eat the most hake?
library(dplyr)
library(tidyr)
library(stats4)
library(ggplot2)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

### Get hake prey lengths from diet data --------------------------------------
path <- "data/diet/CCTD/v4/"

# Read in predator/prey collection and size data & combine
pred <- read.csv(paste0(path, "predator_information_v4.csv"))  # predator data
collection_pred <- read.csv(paste0(path, "collection_information_v4.csv")) %>%
  filter(Collection_ID %in% pred$Collection_ID)
all_pred <- merge(collection_pred, pred, all.y = TRUE)[, c("Month", 
                                                           "Year", 
                                                           "Latitude", 
                                                           "Longitude",
                                                           "Predator_ID",
                                                           "Predator_Com_Name")]


prey_comp <- read.csv(paste0(path, "prey_composition_v4.csv"))
prey_size <-  read.csv(paste0(path, "prey_size_v4.csv")) %>%
  filter(Prey_Comp_ID %in% prey_comp$Prey_Comp_ID)
all_prey <- merge(prey_comp, prey_size, all.y = TRUE)[, c("Prey_Com_Name", 
                                                          "Predator_ID", 
                                                          "Prey_Maturity",
                                                          "Prey_Length1", 
                                                          "Prey_Length_BC_mm")]

pred_prey <- merge(all_pred, all_prey, by = "Predator_ID")

# Filter by hake as prey item, combine length columns into one
hake_preds <-  pred_prey %>% 
  filter(Prey_Com_Name == "Pacific Hake") %>%
  # Make NAs in length columns 0 so they can be added together to combine them
  mutate(Prey_Length1 = replace_na(Prey_Length1, 0),
         Prey_Length_BC_mm = replace_na(Prey_Length_BC_mm, 0)) %>%
  mutate(Prey_Length = Prey_Length1 + Prey_Length_BC_mm)

# Keep only necessary columns and bring NAs back - truly no length recorded
hake_preds <- hake_preds[, c("Year", "Month", "Latitude", "Longitude",
                             "Predator_Com_Name", "Prey_Com_Name",
                             "Prey_Maturity", "Prey_Length")]
hake_preds$Prey_Length[hake_preds$Prey_Length == 0] <- NA


### Convert prey lengths to ages ----------------------------------------------
# Code is further explained in diet_ageing.R
estimate_ages <- function() {
  # Read in age data from the survey & filter to only be hake
  survey_ages <- read.csv("data/diet/acoustic_survey.csv") %>% 
    filter(common_name == "Pacific hake")
  
  age_length <- na.omit(cbind.data.frame(survey = survey_ages$survey,
                                         age = survey_ages$age, 
                                         length = survey_ages$length)) %>%
    mutate(Year = substr(survey, 1, 4))
  
  # Run von Bert using ages from survey data
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
  
  # Calculate prey ages and add to hake predator dataframe
  prey_ages <- (-log(1 - (hake_preds$Prey_Length / 10) / Linf) / K) + a0
  max(na.omit(prey_ages))  # check maximum
  min(na.omit(prey_ages))  # check minimum
  prey_ages[prey_ages < 1] <- 1  # replace values < 1 (lower accumulation age) 
  prey_ages <- round(prey_ages, digits = 0)  # round to whole number
  
  hake_preds <- cbind.data.frame(hake_preds, Estimated_Age = prey_ages) 
  
  hake_preds$Estimated_Age[hake_preds$Prey_Maturity == "Immature"] <- 1
  hake_preds <- hake_preds[!is.na(hake_preds$Estimated_Age), ]
  
  # Read in FEAT data and compare fit 
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
  
  # Plot fit 
  all_ages <- as.data.frame(rbind(cbind(age = hake_preds$Estimated_Age,
                                        length = (hake_preds$Prey_Length / 10),  # prey are in mm
                                        data = rep("CCTD")),
                                  cbind(age = FEAT_hake$predator_age,
                                        length = FEAT_hake$predator_length_cm,
                                        data = rep("FEAT")),
                                  cbind(age = FEAT_hake$prey_ages,
                                        length = FEAT_hake$measure_value,
                                        data = rep("FEAT")))) 
  all_ages$length <- as.numeric(all_ages$length)
  all_ages <- na.omit(all_ages)
  all_ages$data <- factor(all_ages$data, levels = c("CCTD", "FEAT"))
  all_ages$age <- factor(all_ages$age, levels = 1:20)
  
  growth_curve <- ggplot(all_ages, aes(x = age, y = length, color = data)) +
    geom_boxplot(position = "dodge") +
    # scale_alpha_discrete(range = c(0.2, 1)) +
    scale_color_viridis(discrete = TRUE, direction = -1, end = 0.6) +
    xlab("Age") + ylab("Length (cm)") + labs(color = "Source")

  return(list(df = hake_preds, curve = growth_curve))
}
aged <- estimate_ages()$df
estimate_ages()$curve


### Summarize ages consumed ---------------------------------------------------
pred_median_age <- aged %>%
  group_by(Predator_Com_Name) %>%
  summarize(median = median(Estimated_Age),
            mean = mean(Estimated_Age))

consumed_age <- aged %>%
  group_by(Estimated_Age) %>%
  summarize(n = n())
