library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

remove.packages("Rceattle")
devtools::install_github("grantdadams/Rceattle", ref = "dev")

# Read in CEATTLE data from the excel file
hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230324.xlsx")

# Run and fit the CEATTLE model -----------------------------------------------
run_CEATTLE <- function(data, M1, prior, init, msm) {
  data$est_M1 <- M1  # Set M1 to fixed (0), estimated (1), age-varying estimate (3)
  # data$endyr <- 2019
  run <- Rceattle::fit_mod(data_list = data,
                           inits = init,
                           file = NULL, # Don't save
                           msmMode = msm, # Single-species mode - no predation mortality
                           M1Fun = Rceattle::build_M1(M1_model = M1, 
                                                      updateM1 = TRUE,
                                                      M1_use_prior = prior,
                                                      M1_prior_mean = 0.2,
                                                      M1_prior_sd = .1),
                           # proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
                           estimateMode = 0,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
                           HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                                                     FsprLimit = 0.4, # F40%
                                                     Ptarget = 0.4, # Target is 40% B0
                                                     Plimit = 0.1, # No fishing when SB<SB10
                                                     Pstar = 0.45,
                                                     Sigma = 0.5),
                          #  HCR = Rceattle::build_hcr(HCR = 0),
                           phase = "default",
                           initMode = 1)
  
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  return(list(run, fit))
}

# Run in single-species mode
nodiet_fixed <- run_CEATTLE(data = hake_intrasp, M1 = 0, prior = FALSE, init = NULL, msm = 0)
nodiet_fixed[[2]]  # check convergence

nodiet_est <- run_CEATTLE(data = hake_intrasp, M1 = 1, prior = FALSE, init = NULL, msm = 0)
nodiet_est[[2]]  # check convergence

nodiet_prior <- run_CEATTLE(data = hake_intrasp, M1 = 1, prior = TRUE, init = NULL, msm = 0)
nodiet_prior[[2]]  # check convergence

Rceattle::plot_biomass(nodiet_fixed[[1]], add_ci = TRUE)
Rceattle::plot_biomass(nodiet_est[[1]], add_ci = TRUE)
Rceattle::plot_biomass(nodiet_prior[[1]], add_ci = TRUE)


load("main_model.Rdata")

inits <- nodiet_fixed[[1]]$estimated_params
inits$rec_pars <- matrix(1, nrow = nodiet_fixed$data_list$nspp, ncol = 2)
inits$rec_pars[,1] <- inits$ln_mean_rec

inits$ln_mean_rec <- NULL

##########################################
ss_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                            file = NULL,
                            inits = inits, # Initial parameters = 0
                            estimateMode = 4, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1,
                            phase = "default")


nodiet <- nodiet_est[[1]]
# Run with cannibalism, estimated M1
intrasp <-  run_CEATTLE(data = hake_intrasp, M1 = 1, prior = FALSE, init = nodiet$estimated_params, msm = 1)
intrasp[[2]]  # check convergence



# Rceattle diagnostic plots 
# Rceattle::plot_biomass(intrasp[[1]], add_ci = TRUE)
# Rceattle::plot_index(intrasp[[1]])
# Rceattle::plot_catch(intrasp[[1]])
# Rceattle::plot_selectivity(intrasp[[1]])
# Rceattle::plot_mortality(intrasp[[1]], type=3)
# Rceattle::plot_indexresidual(intrasp[[1]])
# Rceattle::plot_logindex(intrasp[[1]])
# Rceattle::plot_recruitment(intrasp[[1]], add_ci = TRUE, incl_proj = TRUE)
# Rceattle::plot_comp(intrasp[[1]])
# Rceattle::plot_srv_comp(intrasp[[1]])
# Rceattle::plot_f(intrasp[[1]])

# Plot population dynamics
start_yr <- intrasp[[1]]$data_list$styr
end_yr <- intrasp[[1]]$data_list$projyr  # end of projection period
hind_end <- intrasp[[1]]$data_list$endyr  # end of hindcast
years <- start_yr:end_yr

plot_models <- function(ms_run, ss_run, assess_yr = as.character(hind_end + 1)) {
  # Plot biomass & recruitment in comparison to no diet & assessment ----------
  ceattle_biomass <- function(run, name) {
    ssb <- (c(run$quantities$biomassSSB) * 2)
    biom <- c(run$quantities$biomass)
    biom_sd <- run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]
    wide <- as.data.frame(cbind(years, ssb, biom))
    colnames(wide) <- c("year", "SSB", "Total Biomass")
    all_biom <- melt(wide, id.vars = "year")
    all_biom2 <- cbind(all_biom, 
                       error = c(run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")], 
                                 run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]), 
                       model = rep(name, length(all_biom$year)))  
    colnames(all_biom2)[2:3] <- c("type", "value")
    return(all_biom2)
  }
  
  biomass <- ceattle_biomass(ms_run, "CEATTLE - cannibalism")
  nodiet_biomass <- ceattle_biomass(ss_run, "CEATTLE - single-species")
  recruitment <- c(ms_run$quantities$R)
  
  # Pull out biomass from stock synthesis & combine, remove pre-1988
  start <- 1966  # start year of SS3 analysis
  ss3_ssb <- cbind(read.table(paste0("data/assessment/", assess_yr, "/ssb.txt"))[-(1:(start_yr-start)), 2:3], 
                   type = rep("SSB"))
  ss3_biomass <- cbind(read.table(paste0("data/assessment/", assess_yr, "/biomass.txt"))[-(1:(start_yr-start)), 2:3], 
                       type = rep("Total Biomass"))
  ss3_biom <- as.data.frame(cbind(year = rep(years, 2), 
                                  rbind(ss3_ssb, ss3_biomass),
                                  model = rep("Assessment")))
  colnames(ss3_biom)[2:3] <- c("value", "error")
  ss3_biom <- ss3_biom[, c(1, 4, 2, 3, 5)]
  
  # Pull out recruitment
  nodiet_R <- c(ss_run$quantities$R)
  ss3_R <- read.table(paste0("data/assessment/", assess_yr, "/recruitment.txt"))[-(1:(start_yr-start)), ]
  
  # Plot biomass & recruitment ------------------------------------------------
  # Put biomass together
  nodiet_biom <- ceattle_biomass(ss_run, "CEATTLE - single-species")
  
  # Combine all biomass sources together
  biom_all <- rbind(biomass, nodiet_biom, ss3_biom)
  biom_all$value <- biom_all$value / 1000000  # to mt
  biom_all$error <- biom_all$error / 1000000  # to mt
  
  # Put recruitment together
  R_wide <- data.frame(year = years, recruitment, nodiet_R)
  colnames(R_wide)[2:3] <- c("CEATTLE - cannibalism", "CEATTLE - single-species")
  R <- melt(R_wide, id.vars = "year")
  # Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS3 is 0)
  ss3_1 <- as.data.frame(cbind(year = ((start_yr+1):end_yr), 
                               variable = rep("Assessment"), 
                               value = ss3_R[1:length((start_yr+1):end_yr), 2],
                               error = ss3_R[1:length((start_yr+1):end_yr), 3]))
  R_all <- rbind(cbind(R, error = c(ms_run$sdrep$sd[which(names(ms_run$sdrep$value) == "R")], 
                                    ss_run$sdrep$sd[which(names(ss_run$sdrep$value) == "R")])), 
                 ss3_1)
  R_all$value <- as.numeric(R_all$value)
  R_all$year <- as.numeric(R_all$year)
  R_all$error <- as.numeric(R_all$error)
  # Reshape to match biomass
  R_new <- cbind(year = R_all$year, 
                 type = rep("Recruitment"), 
                 value = R_all$value / 1000000,  # to millions
                 error = R_all$error / 1000000,  # to millions
                 model = as.character(R_all$variable))
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(biom_all, R_new)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value)
  all_popdy$error <- as.numeric(all_popdy$error)
  all_popdy$model <- factor(all_popdy$model, 
                            levels = c("Assessment", "CEATTLE - single-species", "CEATTLE - cannibalism"))
  all_popdy$type <- factor(all_popdy$type, 
                           labels = c("SSB (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line() +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    geom_vline(xintercept = hind_end, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside") 
  return(popdy_plot)
}

plot_models(intrasp[[1]], nodiet)

# Save plots to specific testing/sensitivity folder
path <- "plots/CEATTLE/cannibalism/Testing/HCR/"
name <- "popdyn_M1prior_HCR6.png"
ggsave(filename=paste0(path, name), 
       plot_models(intrasp[[1]], nodiet[[1]]), 
       width=140, height=150, units="mm", dpi=300)
