# Function to find fit of CEATTLE model
fit <- function(run) {
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  return(fit)
}

# Set up to run Rceattle on main branch ---------------------------------------
remove.packages("Rceattle")
remove.packages("00LOCK-Rceattle")
devtools::install_github("grantdadams/Rceattle")

HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                          FsprLimit = 0.4, # F40%
                          Ptarget = 0.4, # Target is 40% B0
                          Plimit = 0.1, # No fishing when SB<SB10
                          Pstar = 0.45,
                          Sigma = 0.5)

main_data <- Rceattle::read_data(file = "data/hake_intrasp_230427.xlsx")
main_data$est_M1 <- 1
main_run <- Rceattle::fit_mod(data_list = main_data,
                              msmMode = 0, # Single-species mode - no predation mortality
                              estimateMode = 0,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
                              HCR = HCR,
                              phase = "default",
                              projection_uncertainty = TRUE)
fit(main_run)
save(main_run, file = "main_run.Rdata")

###############
## RESTART R ##
###############

# Set up to run Rceattle on dev branch ----------------------------------------
remove.packages("Rceattle")
remove.packages("00LOCK-Rceattle")
devtools::install_github("grantdadams/Rceattle", ref = "dev")

HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                          FsprLimit = 0.4, # F40%
                          Ptarget = 0.4, # Target is 40% B0
                          Plimit = 0.1, # No fishing when SB<SB10
                          Pstar = 0.45,
                          Sigma = 0.5)

dev_data <- Rceattle::read_data(file = "data/hake_intrasp_230427.xlsx")
dev_run <- Rceattle::fit_mod(data_list = dev_data,
                             msmMode = 0, # Single-species mode - no predation mortality
                             estimateMode = 0,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
                             HCR = HCR,
                             phase = "default",
                             projection_uncertainty = TRUE,
                             initMode = 1)

fit(dev_run)

save(dev_run, file = "dev_run.Rdata")


# Compare models across Rceattle verisons -------------------------------------
load("main_run.Rdata")
load("dev_run.Rdata")

# We can plot all runs
mod_list <- list(main_run, dev_run)
mod_names <- c("Main model", "Dev model")

# Plot biomass trajectory
Rceattle::plot_biomass(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, incl_proj = TRUE)


# Can't remember what this is -------------------------------------------------
# inits <- main_run$estimated_params
# inits$rec_pars <- matrix(9, nrow = main_run$data_list$nspp, ncol = 2)
# inits$rec_pars[,1] <- inits$ln_mean_rec
# 
# inits$ln_mean_rec <- NULL

# ##########################################
# ss_run <- Rceattle::fit_mod(data_list = dev_data,
#                             file = NULL,
#                             inits = inits, # Initial parameters = 0
#                             estimateMode = 3, # Estimate
#                             random_rec = FALSE, # No random recruitment
#                             msmMode = 0, # Single species mode
#                             verbose = 1,
#                             phase = "default",
#                             M1Fun = Rceattle::build_M1(M1_model = 1, 
#                                                        updateM1 = FALSE),
#                             initMode = 1)
# 
# ss_run$quantities$jnll_comp
