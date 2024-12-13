#' Script for running the CEATTLE model (with important inputs explicitly listed)
#' for testing the model through running a single- and a multispecies model and
#' plotting the comparison. This script includes a streamlined version of the 
#' population dynamics plots included in ceattle_cannibalism.R and the predation
#' mortality for the multispecies model, as requested by the hake stock 
#' assessors.

# Set up ----------------------------------------------------------------------
# Installation options for Rceattle if testing different branches 
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
library(viridis)
library(here)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Read in CEATTLE data from the excel file
hake_data <- read_data(file = "data/hake_yr24_241126.xlsx")
# start_yr <- new_ms$model$data_list$styr
start_yr <- 1966

# Run and fit the CEATTLE model -----------------------------------------------
run_CEATTLE <- function(data, M1, prior, init, initMode, msm, estMode) {
  data$est_M1 <- M1  
  data$styr <- start_yr
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
                 initMode = initMode,
                 projection_uncertainty = TRUE,
                 random_rec = FALSE,
                 suit_styr = 1993,
                 suit_endyr = 2019) 
  
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
                      estMode = 1,
                      initMode = 1)
new_ms$fit  # check convergence
new_ms$model$quantities$M1
# save(new_ms, file = "models/2024/new_ms_Oct25.Rdata")

# Run single-species model
new_ss <- run_CEATTLE(data = hake_data, 
                      M1 = 0, 
                      prior = TRUE, 
                      init = NULL, 
                      msm = 0, 
                      estMode = 1,
                      initMode = 1)
# new_ss$fit  # check convergence
# save(new_ss, file = "models/2024/new_ss_Oct25.Rdata")

# plot_biomass(Rceattle = list(new_ss$model, new_ms$model), model_names = c("SS", "MS"), add_ci = TRUE)
# plot_ssb(Rceattle = list(new_ss$model, new_ms$model), model_names = c("SS", "MS"), add_ci = TRUE)

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
  
# Plot multispecies vs. single-species ----------------------------------------
end_yr <- 2027
years <- start_yr:end_yr
all_yrs <- start_yr:new_ms$model$data_list$projyr
hind_end <- 2023
plot_start <- 1966
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
  biomass <- ceattle_biomass(ms_run, "cannibalism")
  nodiet_biomass <- ceattle_biomass(ss_run, "single-species")
  
  # Pull out & format recruitment
  nodiet_R <- c(ss_run$quantities$R[, 1:length(start_yr:end_yr)])
  recruitment <- c(ms_run$quantities$R[, 1:length(start_yr:end_yr)])

  # Put biomass together
  nodiet_biom <- ceattle_biomass(ss_run, "single-species")

  # Combine all biomass sources together
  biom_all <- rbind(biomass, nodiet_biom)
  biom_all$value <- biom_all$value / 1000000  # to Mt
  biom_all$error <- biom_all$error / 1000000  # to Mt

  # Put recruitment together
  R_wide <- data.frame(year = years, recruitment, nodiet_R)
  colnames(R_wide)[2:3] <- c("cannibalism", "single-species")
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
  
  # Combine biomass, recruitment
  all_popdy <- rbind(biom_all, R_new)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value)
  all_popdy$error <- as.numeric(all_popdy$error)
  all_popdy$model <- factor(all_popdy$model,
                            levels = c("single-species", 
                                       "cannibalism"))
  all_popdy$type <- factor(all_popdy$type,
                           labels = c("Spawning Biomass (Mt)", 
                                      "Total Biomass (Mt)", 
                                      "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  
  # Plot population dynamics --------------------------------------------------
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value)) +
    geom_line(aes(color = model)) +
    geom_ribbon(aes(ymin=min, ymax=max, fill = model), alpha = 0.2, color = NA) + 
    geom_vline(xintercept = 2024, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    # xlim(plot_start, NA) +
    ylab(" ") + xlab("Year") +
    labs(color = "Model", fill = "Model") +
    facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    ggtitle(title)
    
  return(list(popdy = popdy_plot))
}

plots <- plot_models(new_ms$model, new_ss$model)
plots$popdy

ggsave(plots$popdy,
       file = here("M2", "popdy_241126.png"),
       width=150, height=180, units="mm")

plot_selectivity(new_ms$model)

# Pull out mortality ----------------------------------------------------------
mortality <- function(run, type) {
  M1 <- run$quantities$M1[1, 1, 1:max_age]
  if(type == "single-species") {
    return(M1)
  }
  
  if(type == "multi-species") {
    M2 <- extract_byage(run$quantities$M2, "multispecies", "M2")
    total_mortality <- M2 %>%
      mutate(M1_M2 = M2 + rep(M1, length(all_yrs)))
    total_mortality$age <- as.integer(total_mortality$age)
    total_mortality$year <- as.integer(as.character(total_mortality$year))
    total_mortality <- total_mortality[total_mortality$year <= end_yr,]
    
    mortality_plot <- ggplot(total_mortality, aes(y = age, x = year, zmin = 0, zmax = 1.6)) + 
      geom_tile(aes(fill = M1_M2)) +
      scale_y_continuous(expand = c(0, 0), breaks=c(1, 5, 10, 15, 20)) + 
      scale_x_continuous(expand = c(0, 0)) + 
      scale_fill_viridis(name = "M1 + M2", limits = c(0, 1.5), breaks = c(0.21, 1, 1.5)) +
      geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
      coord_equal() +
      ylab("Age") + xlab("Year") +
      theme(panel.border = element_rect(colour = NA, fill = NA))
    mortality_plot
    
    # Min, max, mean natural mortality by age, for ages 1-5
    M_byage <- total_mortality %>%
      filter(age < 6) %>%
      group_by(age) %>%
      summarize(min = min(M1_M2), max = max(M1_M2), mean = mean(M1_M2))
    
    return(list(M2 = M2, plot = mortality_plot, M_byage = M_byage, M1 = M1, total_M = total_mortality))
  }
}

ms_mort <- mortality(new_ms$model, type = "multi-species")
ms_mort$plot  # plot
ms_mort$M_byage
ms_M1 <- ms_mort$M1
ms_totM <- ms_mort$total_M %>% 
  group_by(age, year) %>%
  summarize(M1_M2 = sum(M1_M2))

ms_M2 <- ms_mort$M2[, -4]
write.csv(ms_M2, here("M2", "M2_241126.csv"), row.names = FALSE)

ms_totM_age1 <- ms_totM %>%
  filter(age == 1)

# # Compare M2 between equilibrium & non-equilibrium models
# m2_init1 <- read.csv(here("M2", "M2_241025_init1.csv"))
# m2_init1$model <- "Eq, initMode = 1"
# m2_init2 <- read.csv(here("M2", "M2_241025_init2.csv"))
# m2_init2$model <- "Non, initMode = 2"
# 
# m2_diff <- rbind.data.frame(m2_init1, m2_init2) %>%
#   filter(year <= 2024) %>%
#   filter(age == 1) %>%
#   ggplot(., aes(x = year, y = M2, color = model)) +
#   geom_line() +
#   ylab("Age 1 M2")
# m2_diff
