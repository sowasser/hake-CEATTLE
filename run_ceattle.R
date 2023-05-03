#' Script for running the CEATTLE model under various configurations for single
#' or multispecies; fixed, estimated, prior M1. These models can take a little 
#' while to run (especially with convergence issues), so output is saved for 
#' use in ceattle_cannibalism.R

library(Rceattle)
# Reinstall Rceattle if needed
# remove.packages("Rceattle")
# pak::pkg_install("grantdadams/Rceattle")
# remotes::install_github("grantdadams/Rceattle")

# Read in CEATTLE data from the excel file
hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230427.xlsx")

### Run and fit the CEATTLE model ---------------------------------------------
run_CEATTLE <- function(data, M1, prior, init, msm, M_phase) {
  data$est_M1 <- M1  
  # data$endyr <- 2019
  run <- fit_mod(data_list = data,
                 inits = init,
                 file = NULL, # Don't save
                 msmMode = msm, # Single-species mode - no predation mortality
                 # M1Fun = Rceattle::build_M1(M1_model = M1, 
                 #                            updateM1 = TRUE,
                 #                            M1_use_prior = prior,
                 #                            M1_prior_mean = 0.2,
                 #                            M1_prior_sd = .1),
                 # proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
                 estimateMode = 0,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
                 HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                                           FsprLimit = 0.4, # F40%
                                           Ptarget = 0.4, # Target is 40% B0
                                           Plimit = 0.1, # No fishing when SB<SB10
                                           Pstar = 0.45,
                                           Sigma = 0.5),
                 M1Fun = Rceattle::build_M1(M1_model = M1, # Set M1 to fixed (0) or estimated (1)
                                            updateM1 = TRUE),
                 phase = "default",
                # # Update phase to help convergence --------------------------
                # phase = list(
                #   dummy = 1,
                #   ln_pop_scalar = 4,
                #   ln_mean_rec = 1,
                #   ln_rec_sigma = 2,
                #   rec_dev = 2,
                #   init_dev = 2,
                #   ln_sex_ratio_sigma = 3,
                #   ln_M1 = M_phase,
                #   ln_mean_F = 1,
                #   ln_Flimit = 3,
                #   ln_Ftarget = 3,
                #   proj_F_prop = 1,
                #   F_dev = 1,
                #   ln_srv_q = 3,
                #   # srv_q_pow = 4,
                #   ln_srv_q_dev = 5,
                #   ln_sigma_srv_q = 4,
                #   ln_sigma_time_varying_srv_q = 4,
                #   sel_coff = 3,
                #   sel_coff_dev = 4,
                #   ln_sel_slp = 3,
                #   sel_inf = 3,
                #   ln_sel_slp_dev = 5,
                #   sel_inf_dev = 5,
                #   ln_sigma_sel = 4,
                #   sel_curve_pen = 4,
                #   ln_sigma_srv_index = 2,
                #   ln_sigma_fsh_catch = 2,
                #   comp_weights = 4,
                #   logH_1 = 6,
                #   logH_1a = 6,
                #   logH_1b = 6,
                #   logH_2 = 6,
                #   logH_3 = 6,
                #   H_4 = 6,
                #   log_gam_a = 5,
                #   log_gam_b = 5,
                #   log_phi = 5),
                # # -----------------------------------------------------------
                # verbose = 1,
                initMode = 1,
                projection_uncertainty = TRUE) 
  
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  
  jnll_summary <- as.data.frame(run$quantities$jnll_comp)
  jnll_summary$sum <- rowSums(run$quantities$jnll_comp)
  
  return(list(model = run, fit = fit, summary = jnll_summary))
}

# Run in single-species mode --------------------------------------------------
ss_fixM1 <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 0, 
                        prior = FALSE, 
                        init = NULL, 
                        msm = 0, 
                        M_phase = 6)
ss_fixM1$fit  # check convergence
save(ss_fixM1, file = "models/ss_fixM1.Rdata")

ss_estM1 <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 1, 
                        prior = FALSE, 
                        init = ss_fixM1[[1]]$estimated_params, 
                        msm = 0, 
                        M_phase = 1)
ss_estM1$fit  # check convergence
ss_estM1$model$quantities$M1
save(ss_estM1, file = "models/ss_estM1.Rdata")

# ss_priorM1 <- run_CEATTLE(data = hake_intrasp, 
#                           M1 = 1, 
#                           prior = TRUE, 
#                           init = NULL, 
#                           msm = 0, 
#                           M_phase = 1)
# ss_priorM1$fit  # check convergence
# ss_priorM1$model$quantities$M1

# Run with cannibalism (multi-species mode) -----------------------------------
ms_estM1 <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 1, 
                        prior = FALSE, 
                        init = ss_estM1$model$estimated_params, 
                        msm = 1, 
                        M_phase = 1)
ms_estM1$fit  # check convergence
# Rceattle diagnostic plots 
Rceattle::plot_biomass(ms_estM1$model, add_ci = TRUE)
Rceattle::plot_recruitment(ms_estM1$model, add_ci = TRUE, incl_proj = TRUE)
save(ms_estM1, file = "models/ms_estM1.Rdata")

ms_fixM1 <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 0, 
                        prior = FALSE, 
                        init = ss_fixM1$model$estimated_params, 
                        msm = 1, 
                        M_phase = 6)
ms_fixM1$fit  # check convergence
save(ms_fixM1, file = "models/ms_fixM1.Rdata")


# For future reference: age-blocked mortality at 1, 3, and 3+ -----------------
# map <- intrasp_run$map
# orig <- intrasp_run$map$mapList$ln_M1
# orig[orig > 2] <- 3
# map$mapList$ln_M1 <- orig
# map$mapFactor$ln_M1 = as.factor(map$mapList$ln_M1)
# 
# intrasp_run2 <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                   inits = NULL,
#                                   file = NULL, # Don't save
#                                   map = map,
#                                   msmMode = 1, # Single-species mode - no predation mortality
#                                   phase = "default")
