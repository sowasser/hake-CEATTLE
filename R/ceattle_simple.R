# Simple script for running the most recent version of the CEATTLE model

library(here)
# devtools::install_github("grantdadams/Rceattle", ref = "dev_srr")
library(Rceattle)

# Read in data
hake_data <- read_data(file = here("data", "hake_yr24_241220.xlsx"))

# Function to collate model output & check convergence
fit_out <- function(run) {
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  return(fit)
}

# Fit the single-species model ------------------------------------------------
ss_run <- fit_mod(data_list = hake_data,
                  inits = NULL,
                  file = NULL, # Don't save
                  msmMode = 0, # Single-species mode - no predation mortality
                  M1Fun = build_M1(M1_model = 0,
                                   updateM1 = TRUE,
                                   M1_use_prior = TRUE,
                                   M_prior = 0.22,
                                   M_prior_sd = .31),
                  # proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
                  estimateMode = 1,  # 0 = Fit the hindcast model and projection with HCR specified via HCR; 1 = hindcast only
                  HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                  DynamicHCR = FALSE,
                                  FsprLimit = 0.4, # F40%
                                  Ptarget = 0.4, # Target is 40% B0
                                  Plimit = 0.1, # No fishing when SB<SB10
                                  Pstar = 0.5,
                                  Sigma = 0.5),
                  phase = TRUE,
                  initMode = 2,  # Population initialization / how to estimate equilibrium age structure
                  projection_uncertainty = TRUE,
                  random_rec = FALSE,
                  suit_styr = 1993,
                  suit_endyr = 2019) 
fit_out(ss_run)  # Check convergence

# Fit the multispecies model --------------------------------------------------
ms_run <- fit_mod(data_list = hake_data,
                  inits = ss_run$model$initial_params,
                  file = NULL, # Don't save
                  msmMode = 1, # Single-species mode - no predation mortality
                  M1Fun = build_M1(M1_model = 1,
                                   updateM1 = TRUE,
                                   M1_use_prior = TRUE,
                                   M_prior = 0.22,
                                   M_prior_sd = .31),
                  # proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
                  estimateMode = 1,  # 0 = Fit the hindcast model and projection with HCR specified via HCR; 1 = hindcast only
                  HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                  DynamicHCR = FALSE,
                                  FsprLimit = 0.4, # F40%
                                  Ptarget = 0.4, # Target is 40% B0
                                  Plimit = 0.1, # No fishing when SB<SB10
                                  Pstar = 0.5,
                                  Sigma = 0.5),
                  phase = TRUE,
                  initMode = 2,  # Population initialization / how to estimate equilibrium age structure
                  projection_uncertainty = TRUE,
                  random_rec = FALSE,
                  suit_styr = 1993,
                  suit_endyr = 2019) 
fit_out(ss_run)  # Check convergence
#' SNW Jan 10: I'm getting a 'convergence warning (8): discontinous likelihood'
#' warning here. If you look at the output, the obj and jnll match and the 
#' gradient is low (enough). I would get this error when there was a slight
#' difference between the obj and the jnll. Not sure why it's happening here,
#' but I don't think it's anything to worry about at the moment.

# Plot models using Rceattle's built-in plotting ------------------------------
plot_biomass(Rceattle = list(ss_run, ms_run), 
             model_names = c("single species", "multispecies"), 
             add_ci = TRUE)
plot_ssb(Rceattle = list(ss_run, ms_run), 
         model_names = c("single species", "multispecies"), 
         add_ci = TRUE)
