# remove.packages("Rceattle")
# remove.packages("00LOCK-Rceattle")
# devtools::install_github("grantdadams/Rceattle", ref = "dev")

# # Local install
# devtools::install_local("~/Desktop/Local/Rceattle")
# library(Rceattle, lib.loc = "~/Desktop/Local/Rceattle")

library(Rceattle)
library(dplyr)
library(reshape2)
library(ggplot2)
library(here)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Read in CEATTLE data from the excel file
hake_data <- read_data(file = "data/hake_yr24_240709.xlsx")

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
                                            M1_prior_mean = 0.22,
                                            M1_prior_sd = .31),
                 # proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
                 estimateMode = estMode,  # 0 = Fit the hindcast model and projection with HCR specified via HCR; 1 = hindcast only
                 HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                                           DynamicHCR = FALSE,
                                           FsprLimit = 0.4, # F40%
                                           Ptarget = 0.4, # Target is 40% B0
                                           Plimit = 0.1, # No fishing when SB<SB10
                                           Pstar = 0.5,
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
                 random_rec = FALSE) 
  
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  
  # # Get table of JNLL components
  # comp <- data.frame(run$quantities$jnll_comp)
  # comp$component <- rownames(comp)
  # rownames(comp) <- NULL
  # comp[nrow(comp) + 1, ] <- c(sum(comp[, 1]), sum(comp[, 2]), "Total NLL")  # add total NLL
  # # Separate comps for fishery & survey (in different columns originally)
  # comp[nrow(comp) + 1, ] <- c(comp[3, 1], 0, "Fishery age composition")
  # comp[nrow(comp) + 1, ] <- c(comp[3, 2], 0, "Survey age composition")
  # comp$NLL <- as.numeric(comp$Sp.Srv.Fsh_1) + as.numeric(comp$Sp.Srv.Fsh_2)  # combine species together
  # comp <- comp[, c("component", "NLL")]
  # comp <- comp %>% filter(NLL != 0)  # remove components w/ no likelihood
  # comp$NLL <- as.numeric(comp$NLL)
  
  return(list(model = run, fit = fit))
}

# Run new model (with cannibalism)
new_ms <- run_CEATTLE(data = hake_data, 
                      M1 = 1, 
                      prior = TRUE, 
                      init = NULL, 
                      msm = 1, 
                      estMode = 0)
new_ms$fit  # check convergence
# new_ms$model$quantities$M1
# save(new_ms, file = "models/2024/new_ms_Jul8.Rdata")

# Run single-species model
new_ss <- run_CEATTLE(data = hake_data, 
                      M1 = 0, 
                      prior = TRUE, 
                      init = NULL, 
                      msm = 0, 
                      estMode = 1)
# new_ss$fit  # check convergence
# save(new_ss, file = "models/2024/new_ss_Jul8.Rdata")

# # Compare to base (publication) model
# load("models/ms_priorM1.Rdata")
# 
# plot_biomass(Rceattle = list(new_ms$model, ms_priorM1$model),
#              model_names = c("New Model", "Base Model"),
#              incl_proj = TRUE,
#              add_ci = TRUE)
# plot_ssb(Rceattle = list(new_ms$model, ms_priorM1$model),
#          model_names = c("New Model", "Base Model"),
#          incl_proj = TRUE,
#          add_ci = TRUE)
  
### Plot multi-species vs. single-species vs. assessment ----------------------
start_yr <- new_ms$model$data_list$styr
end_yr <- 2024
years <- start_yr:end_yr
all_yrs <- new_ms$model$data_list$styr:new_ms$model$data_list$projyr
hind_end <- 2023
max_age <- new_ms$model$data_list$nages

# Helper function for extracting -by-age data from CEATTLE
extract_byage <- function(result, name, type) {
  df <- as.data.frame(as.table(result))

  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(all_yrs)
  colnames(df) <- c("age", "year", type)

  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"

  return(df)
}

plot_models <- function(ms_run, ss_run, title = "") {
  # ms_run <- new_ms$model
  # ss_run <- new_ss$model
  # Helper function to get biomass from CEATTLE
  ceattle_biomass <- function(run, name) {
    ssb <- (c(run$quantities$biomassSSB[, 1:length(start_yr:end_yr)]) * 2)
    biom <- c(run$quantities$biomass[, 1:length(start_yr:end_yr)])
    biom_sd <- run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]
    wide <- as.data.frame(cbind(years, ssb, biom))
    colnames(wide) <- c("year", "SSB", "Total Biomass")
    all_biom <- melt(wide, id.vars = "year")
    all_biom2 <- cbind(all_biom,
                       error = c(run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")][1:length(start_yr:end_yr)],
                                 run$sdrep$sd[which(names(run$sdrep$value) == "biomass")][1:length(start_yr:end_yr)]),
                       model = rep(name, length(all_biom$year)))
    colnames(all_biom2)[2:3] <- c("type", "value")
    return(all_biom2)
  }
  # Pull results from CEATTLE models & assessment -----------------------------
  # Pull out biomass from CEATTLE 
  biomass <- ceattle_biomass(ms_run, "CEATTLE - cannibalism")
  nodiet_biomass <- ceattle_biomass(ss_run, "CEATTLE - single-species")
  
  # Get assessment derived quantities
  assess <- read.csv(here("data", "assessment", "2024", "median-population-estimates.csv"))
  
  # Pull out & format biomass
  assess_bio <- assess[, c(1, 2, 5)]
  colnames(assess_bio) <- c("year", "SSB", "Total Biomass")
  assess_bio[, 2] <- assess_bio[, 2] / 1000  # to millions
  assess_bio[, 2] <- assess_bio[, 2] * 2 # both sexes
  assess_bio[, 3] <- assess_bio[, 3] / 1000  # to millions
  assess_bio <- melt(assess_bio, id.vars = "year")
  colnames(assess_bio)[2] <- "type"
  assess_bio$error <- 0
  assess_bio$model <- "Assessment"
  assess_bio <- assess_bio %>% filter(year %in% years)
  
  # Pull out & format recruitment
  nodiet_R <- c(ss_run$quantities$R[, 1:length(start_yr:end_yr)])
  recruitment <- c(ms_run$quantities$R[, 1:length(start_yr:end_yr)])
  
  assess_rec <- assess[, c(1, 6)]
  assess_rec$Year <- assess_rec$Year + 1
  assess_rec <- assess_rec %>% filter(Year %in% years)
  assess_rec <- cbind(year = assess_rec$Year,
                      type = rep("Recruitment"),
                      value = assess_rec[, 2] / 1000,  
                      error = 0,
                      model = "Assessment")

  # Put biomass together
  nodiet_biom <- ceattle_biomass(ss_run, "CEATTLE - single-species")

  # Combine all biomass sources together
  biom_all <- rbind(biomass, nodiet_biom)
  biom_all$value <- biom_all$value / 1000000  # to Mt
  biom_all$error <- biom_all$error / 1000000  # to Mt
  biom_all <- rbind(biom_all, assess_bio)

  # Put recruitment together
  R_wide <- data.frame(year = years, recruitment, nodiet_R)
  colnames(R_wide)[2:3] <- c("CEATTLE - cannibalism", "CEATTLE - single-species")
  R <- melt(R_wide, id.vars = "year")
  R_all <- rbind(cbind(R, error = c(ms_run$sdrep$sd[which(names(ms_run$sdrep$value) == "R")][1:length(start_yr:end_yr)],
                                    ss_run$sdrep$sd[which(names(ss_run$sdrep$value) == "R")][1:length(start_yr:end_yr)])))
  R_all$value <- as.numeric(R_all$value)
  R_all$year <- as.numeric(R_all$year)
  R_all$error <- as.numeric(R_all$error)
  # Reshape to match biomass
  R_new <- cbind(year = R_all$year,
                 type = rep("Recruitment"),
                 value = R_all$value / 1000000,  # to millions
                 error = R_all$error / 1000000,  # to millions
                 model = as.character(R_all$variable))
  R_new <- rbind(R_new, assess_rec)
  
  # Combine relative SSB together
  assess_relSSB <- assess[, c(1, 3)] %>% filter(Year %in% years)
  relSSB <- rbind.data.frame(cbind.data.frame(depletion = t(ss_run$quantities$depletionSSB)[1:length(years)],
                                              year = years, 
                                              model = "CEATTLE - single-species"),
                             cbind.data.frame(depletion = t(ms_run$quantities$depletionSSB)[1:length(years)],
                                              year = years, 
                                              model = "CEATTLE - cannibalism"),
                             cbind.data.frame(depletion = assess_relSSB[, 2] / 100, # to decimal
                                              year = years,
                                              model = "Assessment"))
  
  # Set columns to match popdy dataframe
  relSSB$error <- 0
  relSSB$type <- "Relative SB"
  relSSB <- relSSB[, c(2, 5, 1, 4, 3)]
  colnames(relSSB)[3] <- "value"

  # Combine biomass, recruitment, relative SB 
  all_popdy <- rbind(biom_all, R_new, relSSB)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value)
  all_popdy$error <- as.numeric(all_popdy$error)
  all_popdy$model <- factor(all_popdy$model,
                            levels = c("Assessment", 
                                       "CEATTLE - single-species", 
                                       "CEATTLE - cannibalism"))
  all_popdy$type <- factor(all_popdy$type,
                           labels = c("Spawning Biomass (Mt)", 
                                      "Total Biomass (Mt)", 
                                      "Recruitment (millions)", 
                                      "Relative Spawning Biomass"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  
  # Plot population dynamics --------------------------------------------------
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value)) +
    geom_hline(data = all_popdy %>% filter(type == "Relative Spawning Biomass"),
               aes(yintercept = 1), col = "gray") +
    geom_hline(data = all_popdy %>% filter(type == "Relative Spawning Biomass"),
               aes(yintercept = 0.4), col = "gray") +
    geom_hline(data = all_popdy %>% filter(type == "Relative Spawning Biomass"),
               aes(yintercept = 0.1), col = "gray") +
    geom_hline(data = all_popdy %>% filter(type == "Relative Spawning Biomass"),
               aes(yintercept = 1.3), col = "white") +  # hack to control y-axis limits for relSSB facet
    geom_line(aes(color = model)) +
    geom_ribbon(aes(ymin=min, ymax=max, fill = model), alpha = 0.2, color = NA) + 
    geom_vline(xintercept = 2024, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    ylab(" ") + xlab("Year") +
    labs(color = "Model", fill = "Model") +
    facet_wrap(~type, ncol = 2, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    ggtitle(title)
  
  
  # Difference between single-species and assessment models -------------------
  # Read in 2020 model runs and assessment results
  load("models/ms_priorM1.rdata")
  load("models/ss_priorM1.rdata")
  assess_yr <- 2020
  biomass_2020 <- ceattle_biomass(ms_priorM1$model, "CEATTLE - cannibalism")
  nodiet_biomass_2020 <- ceattle_biomass(ss_priorM1$model, "CEATTLE - single-species")
  recruitment_2020 <- c(ms_priorM1$model$quantities$R[, 1:length(start_yr:2022)])
  
  start <- 1966  # start year of SS3 analysis
  ss3_ssb_2020 <- cbind(read.table(paste0("data/assessment/", assess_yr, "/ssb.txt"))[-(1:(start_yr-start)), 2:3], 
                        type = rep("SSB"))
  ss3_biomass_2020 <- cbind(read.table(paste0("data/assessment/", assess_yr, "/biomass.txt"))[-(1:(start_yr-start)), 2:3], 
                            type = rep("Total Biomass"))
  ss3_biom_2020 <- as.data.frame(cbind(year = rep(start_yr:2022, 2), 
                                       rbind(ss3_ssb_2020, ss3_biomass_2020),
                                       model = "2020 Assessment"))
  colnames(ss3_biom_2020)[2:3] <- c("value", "error")
  ss3_biom_2020 <- ss3_biom_2020[, c(1, 4, 2, 3, 5)]
  
  # Pull out recruitment
  nodiet_R_2020 <- c(ss_priorM1$model$quantities$R[, 1:length(start_yr:2022)])
  
  # Get recruitment from SS3 files, offset by 1 year
  ss3_R_2020 <- read.table(paste0("data/assessment/", assess_yr, "/recruitment.txt"))[-(1:((start_yr-1)-start)), ]
  ss3_R_2020 <- ss3_R_2020[-nrow(ss3_R_2020), ]
  
  # Combine all biomass sources together
  biom_all_2020 <- rbind(biomass_2020, nodiet_biomass_2020, ss3_biom_2020)
  biom_all_2020$value <- biom_all_2020$value / 1000000  # to Mt
  biom_all_2020$error <- biom_all_2020$error / 1000000  # to Mt
  
  # Put recruitment together
  R_wide_2020 <- data.frame(year = start_yr:2022, recruitment_2020, nodiet_R_2020)
  colnames(R_wide_2020)[2:3] <- c("CEATTLE - cannibalism", "CEATTLE - single-species")
  R_2020 <- melt(R_wide_2020, id.vars = "year")
  # Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS3 is 0)
  ss3_R_2020 <- as.data.frame(cbind(year = start_yr:2022, 
                                    variable = "2020 Assessment", 
                                    value = ss3_R_2020[, 2],
                                    error = ss3_R_2020[, 3]))
  R_all_2020 <- rbind(cbind(R_2020, error = c(ms_priorM1$model$sdrep$sd[which(names(ms_priorM1$model$sdrep$value) == "R")][1:length(start_yr:2022)], 
                                              ss_priorM1$model$sdrep$sd[which(names(ss_priorM1$model$sdrep$value) == "R")][1:length(start_yr:2022)])), 
                 ss3_R_2020)
  R_all_2020$value <- as.numeric(R_all_2020$value)
  R_all_2020$year <- as.numeric(R_all_2020$year)
  R_all_2020$error <- as.numeric(R_all_2020$error)
  # Reshape to match biomass
  R_new_2020 <- cbind(year = R_all_2020$year, 
                      type = rep("Recruitment"), 
                      value = R_all_2020$value / 1000000,  # to millions
                      error = R_all_2020$error / 1000000,  # to millions
                      model = as.character(R_all_2020$variable))
  
  # Combine biomass & recruitment and plot
  all_popdy_2020 <- rbind(biom_all_2020, R_new_2020)
  all_popdy_2020$year <- as.numeric(all_popdy_2020$year)
  all_popdy_2020$value <- as.numeric(all_popdy_2020$value)
  all_popdy_2020$error <- as.numeric(all_popdy_2020$error)
  all_popdy_2020 <- all_popdy_2020 %>% filter(year <= 2022)  # make sure everything only goes to 2022
  
  # Calculate difference for each assessment period
  get_diff <- function(df, stat, assess, assess_yr, new_type) {
    df_assess <- df %>% filter(type == stat & model == assess)
    df_ceattle <- df %>% filter(type == stat & model == "CEATTLE - single-species")
    diff <- cbind.data.frame(year = df$year,
                             type = new_type,
                             value = df_assess$value - df_ceattle$value,
                             Assessment = rep(assess_yr))
    return(diff)
  }
  
  # Combine and plot difference -----------------------------------------------
  diff <- rbind.data.frame(get_diff(all_popdy_2020, "Total Biomass", "2020 Assessment", factor(2020), "Total Biomass"),
                           get_diff(all_popdy, "Total Biomass (Mt)", "Assessment", factor(2024), "Total Biomass"),
                           get_diff(all_popdy_2020, "SSB", "2020 Assessment", factor(2020), "Spawning Biomass"),
                           get_diff(all_popdy, "Spawning Biomass (Mt)", "Assessment", factor(2024), "Spawning Biomass"),
                           get_diff(all_popdy_2020, "Recruitment", "2020 Assessment", factor(2020), "Recruitment"),
                           get_diff(all_popdy, "Recruitment (millions)", "Assessment", factor(2024), "Recruitment")) %>%
    mutate(type = factor(type, levels = c("Spawning Biomass", "Total Biomass", "Recruitment"))) %>%
    ggplot(., aes(x = year, y = value, color = Assessment)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "gray") +
      geom_line() +
      facet_wrap(~type, scales = "free_y") +
      ylab("Assessment - single-species model")
    
  return(list(popdy = popdy_plot, diff = diff))
}

plots <- plot_models(new_ms$model, new_ss$model)
plots$popdy
plots$diff

ggsave(plots$popdy,
       file = here("M2", "2024 assessment", "age1-index.png"),
       width=270, height=150, units="mm")

ggsave(plots$diff,
       file = here("M2", "2024 assessment", "assess_diff.png"),
       width = 270, height = 90, units = "mm")

plot_selectivity(new_ms$model)

# # set up r4ss output (kinda working) ------------------------------------------
# Dirplot <- here("data", "assessment", "2024", "mcmc", "sso")
# 
# replist <- r4ss::SS_output(
#   dir = Dirplot,
#   verbose = TRUE,
#   printstats = FALSE,
#   covar = FALSE
# )
# 
# # plots the results (store in the 'plots' sub-directory)
# r4ss::SS_plots(replist,
#          dir = Dirplot,
#          plot = 2:26,
#          printfolder = 'plots'
# )
