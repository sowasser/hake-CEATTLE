#' Script for running the CEATTLE model under various configurations for single
#' or multispecies; fixed, estimated, prior M1. These models can take a little 
#' while to run (especially with convergence issues), so output is saved for 
#' use in ceattle_cannibalism.R

# Reinstall Rceattle if needed
# remove.packages("Rceattle")
# remove.packages("00LOCK-Rceattle")
# pak::pkg_install("grantdadams/Rceattle")
# remotes::install_github("grantdadams/Rceattle")
# devtools::install_github("grantdadams/Rceattle", ref = "dev")
library(Rceattle)
library(dplyr)
library(scales)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Read in CEATTLE data from the excel file
hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230808_a20.xlsx")

### Run and fit the CEATTLE model ---------------------------------------------
run_CEATTLE <- function(data, M1, prior, init, msm, estMode) {
  data$est_M1 <- M1  
  # data$endyr <- 2019
  run <- fit_mod(data_list = data,
                 inits = init,
                 file = NULL, # Don't save
                 msmMode = msm, # Single-species mode - no predation mortality
                 M1Fun = Rceattle::build_M1(M1_model = M1,
                                            updateM1 = TRUE,
                                            M1_use_prior = prior,
                                            M1_prior_mean = 0.2,
                                            M1_prior_sd = .1),
                 # proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
                 estimateMode = estMode,  # 0 = Fit the hindcast model and projection with HCR specified via HCR; 1 = hindcast only
                 HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                                           FsprLimit = 0.4, # F40%
                                           Ptarget = 0.4, # Target is 40% B0
                                           Plimit = 0.1, # No fishing when SB<SB10
                                           Pstar = 0.45,
                                           Sigma = 0.5),
                 phase = "default",
                 # Update phase to help convergence ---------------------------
                 # phase = list(
                 #   dummy = 1,
                 #   ln_pop_scalar = 4, # Scalar for input numbers-at-age
                 #   rec_pars = 1, # Stock-recruit parameters or log(mean rec) if no stock-recruit relationship
                 #   ln_rec_sigma = 2, # Variance for annual recruitment deviats
                 #   rec_dev = 2, # Annual recruitment deviats
                 #   init_dev = 2, # Age specific initial age-structure deviates or parameters
                 #   ln_sex_ratio_sigma = 3, # Variance of sex ratio (usually fixed)
                 #   ln_M1 = 7, #  Estimated natural or residual mortality
                 #   ln_mean_F = 1, # Mean fleet-specific fishing mortality
                 #   ln_Flimit = 3, # Estimated F limit
                 #   ln_Ftarget = 3, # Estimated F target
                 #   ln_Finit = 3, # Estimated fishing mortality for non-equilibrium initial age-structure
                 #   proj_F_prop = 1, # Fixed fleet-specific proportion of Flimit and Ftarget apportioned within each species
                 #   F_dev = 1, # Annual fleet specific fishing mortality deviates
                 #   ln_srv_q = 3, # Survey catchability
                 #   ln_srv_q_dev = 5, # Annual survey catchability deviates (if time-varying)
                 #   ln_sigma_srv_q = 4, # Prior SD for survey catchability deviates
                 #   ln_sigma_time_varying_srv_q = 4, # SD for annual survey catchability deviates (if time-varying)
                 #   sel_coff = 3, # Non-parametric selectivity coefficients
                 #   sel_coff_dev = 4, # Annual deviates for non-parametric selectivity coefficients
                 #   ln_sel_slp = 3, # Slope parameters for logistic forms of selectivity
                 #   sel_inf = 3, # Asymptote parameters for logistic forms of selectivity
                 #   ln_sel_slp_dev = 5, # Annual deviates for slope parameters for logistic forms of selectivity (if time-varying)
                 #   sel_inf_dev = 5, # Annual deviates for asymptote parameters for logistic forms of selectivity (if time-varying)
                 #   ln_sigma_sel = 4, # SD for annual selectivity deviates (if time-varying)
                 #   sel_curve_pen = 4, # Penalty for non-parametric selectivity
                 #   ln_sigma_srv_index = 2, # Log SD for survey lognormal index likelihood (usually input)
                 #   ln_sigma_fsh_catch = 2, # Log SD for lognormal catch likelihood (usually input)
                 #   comp_weights = 4, # Weights for multinomial comp likelihood
                 #   logH_1 = 6,  # Functional form parameter (not used in MSVPA functional form)
                 #   logH_1a = 6, # Functional form parameter (not used in MSVPA functional form)
                 #   logH_1b = 6, # Functional form parameter (not used in MSVPA functional form)
                 #   logH_2 = 6, # Functional form parameter (not used in MSVPA functional form)
                 #   logH_3 = 6, # Functional form parameter (not used in MSVPA functional form)
                 #   H_4 = 6, # Functional form parameter (not used in MSVPA functional form)
                 #   log_gam_a = 5, # Suitability parameter (not used in MSVPA style)
                 #   log_gam_b = 5, # Suitability parameter (not used in MSVPA style)
                 #   log_phi = 5 # Suitability parameter (not used in MSVPA style)
                 # ),
                 # ------------------------------------------------------------
                 initMode = 1,
                 projection_uncertainty = TRUE,
                 loopnum = 7) 
  
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  
  # Get table of JNLL components
  comp <- data.frame(run$quantities$jnll_comp)
  comp$component <- rownames(comp)
  rownames(comp) <- NULL
  comp[nrow(comp) + 1, ] <- c(sum(comp[, 1]), sum(comp[, 2]), "Total NLL")  # add total NLL
  # Separate comps for fishery & survey (in different columns originally)
  comp[nrow(comp) + 1, ] <- c(comp[3, 1], 0, "Fishery age composition")
  comp[nrow(comp) + 1, ] <- c(comp[3, 2], 0, "Survey age composition")
  comp$NLL <- as.numeric(comp$Sp.Srv.Fsh_1) + as.numeric(comp$Sp.Srv.Fsh_2)  # combine species together
  comp <- comp[, c("component", "NLL")]
  comp <- comp %>% filter(NLL != 0)  # remove components w/ no likelihood
  comp$NLL <- as.numeric(comp$NLL)
  
  return(list(model = run, fit = fit, summary = comp))
}

# Run in single-species mode --------------------------------------------------
ss_fixM1 <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 0, 
                        prior = FALSE, 
                        init = NULL, 
                        msm = 0, 
                        estMode = 0)

ss_fixM1$fit  # check convergence
save(ss_fixM1, file = "models/ss_fixM1.Rdata")

ss_estM1 <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 1, 
                        prior = FALSE, 
                        init = ss_fixM1[[1]]$estimated_params, 
                        msm = 0, 
                        estMode = 0)
ss_estM1$fit  # check convergence
ss_estM1$model$quantities$M1
save(ss_estM1, file = "models/ss_estM1.Rdata")

ss_priorM1 <- run_CEATTLE(data = hake_intrasp,
                          M1 = 1,
                          prior = TRUE,
                          init = NULL,
                          msm = 0,
                          estMode = 0)
ss_priorM1$fit  # check convergence
ss_priorM1$model$quantities$M1
save(ss_priorM1, file = "models/ss_priorM1.Rdata")

# Run with cannibalism (multi-species mode) -----------------------------------
ms_fixM1 <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 0, 
                        prior = FALSE, 
                        init = ss_fixM1$model$estimated_params, 
                        msm = 1, 
                        estMode = 0)
ms_fixM1$fit  # check convergence
save(ms_fixM1, file = "models/ms_fixM1.Rdata")

ms_estM1 <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 1, 
                        prior = FALSE, 
                        init = ms_fixM1$model$estimated_params, 
                        msm = 1, 
                        estMode = 0)
ms_estM1$fit  # check convergence
ms_estM1$model$quantities$M1
# # Rceattle diagnostic plots 
# Rceattle::plot_biomass(ms_estM1$model, add_ci = TRUE)
# Rceattle::plot_recruitment(ms_estM1$model, add_ci = TRUE, incl_proj = TRUE)
save(ms_estM1, file = "models/ms_estM1.Rdata")

ms_priorM1 <- run_CEATTLE(data = hake_intrasp,
                          M1 = 1,
                          prior = TRUE,
                          init = NULL,
                          msm = 1,
                          estMode = 0)
ms_priorM1$fit  # check convergence
ms_priorM1$model$quantities$M1
save(ms_priorM1, file = "models/ms_priorM1.Rdata")

data_noproj <- hake_intrasp
data_noproj$projyr <- 2019
ms_noproj <- run_CEATTLE(data = data_noproj, 
                         M1 = 1, 
                         prior = FALSE, 
                         init = NULL, 
                         msm = 1, 
                         estMode = 1)
ms_noproj$fit  # check convergence
save(ms_noproj, file = "models/ms_noproj.Rdata")
ms_priorM1_noproj <- run_CEATTLE(data = data_noproj,
                                 M1 = 1,
                                 prior = TRUE,
                                 init = NULL,
                                 msm = 1,
                                 estMode = 1)
ms_priorM1_noproj$fit  # check convergence
save(ms_priorM1_noproj, file = "models/ms_priorM1_noproj.Rdata")


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


### Run sensitivities ---------------------------------------------------------
# Time-varying diet -----------------------------------------------------------
# Read in different Dirichlet-corrected datasets
dirichlet_90s <- read.csv("data/diet/Dirichlet/Dirichlet_90s.csv")
data_90s <- hake_intrasp
data_90s$UobsWtAge <- dirichlet_90s
data_90s$styr <- 1980
data_90s$endyr <- 2019
data_90s$projyr <- 2019
run_90s <- run_CEATTLE(data = data_90s,
                       M1 = 1,
                       prior = FALSE,
                       init = NULL,
                       msm = 1,
                       estMode = 1)
run_90s$fit  # check convergence
save(run_90s, file = "models/sensitivity/time-varying/run_90s.Rdata")
run_90s_prior <- run_CEATTLE(data = data_90s,
                             M1 = 1,
                             prior = TRUE,
                             init = NULL,
                             msm = 1,
                             estMode = 1)
run_90s_prior$fit  # check convergence
save(run_90s_prior, file = "models/sensitivity/time-varying/run_90s_prior.Rdata")

dirichlet_recent <- read.csv("data/diet/Dirichlet/Dirichlet_recent.csv")
data_recent <- hake_intrasp
data_recent$UobsWtAge <- dirichlet_recent
data_recent$styr <- 1980
data_recent$endyr <- 2019
data_recent$projyr <- 2019
run_recent <- run_CEATTLE(data = data_recent,
                          M1 = 1,
                          prior = FALSE,
                          init = NULL,
                          msm = 1,
                          estMode = 1)
run_recent$fit  # check convergence
save(run_recent, file = "models/sensitivity/time-varying/run_recent.Rdata")
run_recent_prior <- run_CEATTLE(data = data_recent,
                                M1 = 1,
                                prior = TRUE,
                                init = NULL,
                                msm = 1,
                                estMode = 1)
run_recent_prior$fit  # check convergence
save(run_recent_prior, file = "models/sensitivity/time-varying/run_recent_prior.Rdata")


# Variation in diet proportion ------------------------------------------------
# Pull out data from base intrasp run
wts <- hake_intrasp$UobsWtAge %>% 
  group_by(Pred_age, Prey_age) %>%
  summarize(wt_prop = mean(Stomach_proportion_by_weight))

wt05 <- rescale_max(wts$wt_prop, to = c(0, 0.005))
wt10 <- rescale_max(wts$wt_prop, to = c(0, 0.1))
wt50 <- rescale_max(wts$wt_prop, to = c(0, 0.5))
wt75 <- rescale_max(wts$wt_prop, to = c(0, 0.75))

# # Look at diet proportions
# prop <- as.data.frame(cbind(wts, wt05 = wt05, wt10 = wt10, wt50 = wt50, wt75 = wt75))
# colnames(prop)[3] <- c("observed data")
# prop_all <- reshape2::melt(prop, id.vars = c("Pred_age", "Prey_age"))
# 
# stomach_props <- ggplot(prop_all, aes(x=Prey_age, y=value, fill=variable)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
#   ylab("stomach proportion") + xlab("prey age") +
#   facet_wrap(~Pred_age, ncol = 3)
# stomach_props
# 
# ggsave(filename = "plots/CEATTLE/cannibalism/Testing/sensitivity_prop.png", stomach_props,
#        width=140, height=150, units="mm", dpi=300)

data05 <- hake_intrasp
data05$UobsWtAge$Stomach_proportion_by_weight <- wt05
run_wt05 <- run_CEATTLE(data = data05,
                        M1 = 1,
                        prior = FALSE,
                        init = NULL,
                        msm = 1,
                        estMode = 0)
run_wt05$fit
save(run_wt05, file = "models/sensitivity/diet/run_wt05.Rdata")
run_wt05_fix <- run_CEATTLE(data = data05,
                            M1 = 0,
                            prior = FALSE,
                            init = NULL,
                            msm = 1,
                            estMode = 0)
run_wt05_fix$fit
save(run_wt05_fix, file = "models/sensitivity/diet/run_wt05_fix.Rdata")
run_wt05_prior <- run_CEATTLE(data = data05,
                              M1 = 1,
                              prior = TRUE,
                              init = NULL,
                              msm = 1,
                              estMode = 0)
run_wt05_prior$fit
save(run_wt05_prior, file = "models/sensitivity/diet/run_wt05_prior.Rdata")
data05$projyr <- 2019
run_wt05_noproj <- run_CEATTLE(data = data05,
                               M1 = 1,
                               prior = FALSE,
                               init = NULL,
                               msm = 1,
                               estMode = 1)
run_wt05_noproj$fit
save(run_wt05_noproj, file = "models/sensitivity/diet/run_wt05_noproj.Rdata")

data10 <- hake_intrasp
data10$UobsWtAge$Stomach_proportion_by_weight <- wt10
run_wt10 <- run_CEATTLE(data = data10,
                        M1 = 1,
                        prior = FALSE,
                        init = NULL,
                        msm = 1,
                        estMode = 0)
run_wt10$fit
save(run_wt10, file = "models/sensitivity/diet/run_wt10.Rdata")
run_wt10_fix <- run_CEATTLE(data = data10,
                            M1 = 0,
                            prior = FALSE,
                            init = NULL,
                            msm = 1,
                            estMode = 0)
run_wt10_fix$fit
save(run_wt10_fix, file = "models/sensitivity/diet/run_wt10_fix.Rdata")
run_wt10_prior <- run_CEATTLE(data = data10,
                              M1 = 1,
                              prior = TRUE,
                              init = NULL,
                              msm = 1,
                              estMode = 0)
run_wt10_prior$fit
save(run_wt10_prior, file = "models/sensitivity/diet/run_wt10_prior.Rdata")
data10$projyr <- 2019
run_wt10_noproj <- run_CEATTLE(data = data10,
                               M1 = 1,
                               prior = FALSE,
                               init = NULL,
                               msm = 1,
                               estMode = 1)
run_wt10_noproj$fit
save(run_wt10_noproj, file = "models/sensitivity/diet/run_wt10_noproj.Rdata")

data50 <- hake_intrasp
data50$UobsWtAge$Stomach_proportion_by_weight <- wt50
run_wt50 <- run_CEATTLE(data = data50,
                        M1 = 1,
                        prior = FALSE,
                        init = NULL,
                        msm = 1,
                        estMode = 0)
run_wt50$fit
save(run_wt50, file = "models/sensitivity/diet/run_wt50.Rdata")
run_wt50_fix <- run_CEATTLE(data = data50,
                            M1 = 0,
                            prior = FALSE,
                            init = NULL,
                            msm = 1,
                            estMode = 0)
run_wt50_fix$fit
save(run_wt50_fix, file = "models/sensitivity/diet/run_wt50_fix.Rdata")
run_wt50_prior <- run_CEATTLE(data = data50,
                              M1 = 1,
                              prior = TRUE,
                              init = NULL,
                              msm = 1,
                              estMode = 0)
run_wt50_prior$fit
save(run_wt50_prior, file = "models/sensitivity/diet/run_wt50_prior.Rdata")
data50$projyr <- 2019
run_wt50_noproj <- run_CEATTLE(data = data50,
                               M1 = 1,
                               prior = FALSE,
                               init = NULL,
                               msm = 1,
                               estMode = 1)
run_wt50_noproj$fit
save(run_wt50_noproj, file = "models/sensitivity/diet/run_wt50_noproj.Rdata")

data75 <- hake_intrasp
data75$UobsWtAge$Stomach_proportion_by_weight <- wt75
run_wt75 <- run_CEATTLE(data = data75,
                        M1 = 1,
                        prior = FALSE,
                        init = NULL,
                        msm = 1,
                        estMode = 0)
run_wt75$fit
save(run_wt75, file = "models/sensitivity/diet/run_wt75.Rdata")
run_wt75_fix <- run_CEATTLE(data = data75,
                            M1 = 0,
                            prior = FALSE,
                            init = NULL,
                            msm = 1,
                            estMode = 0)
run_wt75_fix$fit
save(run_wt75_fix, file = "models/sensitivity/diet/run_wt75_fix.Rdata")
run_wt75_prior <- run_CEATTLE(data = data75,
                              M1 = 1,
                              prior = TRUE,
                              init = NULL,
                              msm = 1,
                              estMode = 0)
run_wt75_prior$fit
save(run_wt75_prior, file = "models/sensitivity/diet/run_wt75_prior.Rdata")
data75$projyr <- 2019
run_wt75_noproj <- run_CEATTLE(data = data75,
                               M1 = 1,
                               prior = FALSE,
                               init = NULL,
                               msm = 1,
                               estMode = 1)
run_wt75_noproj$fit
save(run_wt75_noproj, file = "models/sensitivity/diet/run_wt75_noproj.Rdata")
