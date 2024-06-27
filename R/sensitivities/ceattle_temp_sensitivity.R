# Script for testing the cannibalism model's sensitivity to temperature

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230808_a20.xlsx")
load("models/ms_priorM1.Rdata")

# Pull out temperature from base run & create temps for sensitivities
temp_low <- hake_intrasp$env_data$BTempC - 2
temp_high <- hake_intrasp$env_data$BTempC + 2

run_ceattle <- function(temp, df) {
  df$env_data$BTempC <- temp
  ceattle <- Rceattle::fit_mod(data_list = df,
                               inits = ms_priorM1$model$estimated_params, # Initial parameters = 0
                               file = NULL, # Don't save
                               M1Fun = Rceattle::build_M1(M1_model = 1,
                                                          updateM1 = TRUE,
                                                          M1_use_prior = TRUE,
                                                          M1_prior_mean = 0.2,
                                                          M1_prior_sd = .1),
                               HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                                                         FsprLimit = 0.4, # F40%
                                                         Ptarget = 0.4, # Target is 40% B0
                                                         Plimit = 0.1, # No fishing when SB<SB10
                                                         Pstar = 0.45,
                                                         Sigma = 0.5),
                               msmMode = 1, # Multi-species mode
                               phase = "default",
                               initMode = 1,
                               projection_uncertainty = TRUE,
                               loopnum = 7)
  return(ceattle)
}

# Run low-cannibalism models
run_low <- run_ceattle(temp_low, hake_intrasp)
run_high <- run_ceattle(temp_high, hake_intrasp)

Rceattle::plot_biomass(Rceattle = list(ms_priorM1$model, run_low, run_high),
                       model_names = c("cannibalism", "low", "high"))
Rceattle::plot_recruitment(Rceattle = list(ms_priorM1$model, run_low, run_high),
                           model_names = c("cannibalism", "low", "high"))

# # Plot population dynamics across model runs ----------------------------------
# years <- 1988:2019
# 
# test_plot_popdy <- function() {
#   # Pull out SSB & overall biomass from CEATTLE runs
#   ceattle_biomass <- function(run, name) {
#     ssb <- (c(run$quantities$biomassSSB) * 2)
#     biom <- c(run$quantities$biomass)
#     wide <- as.data.frame(cbind(years, ssb, biom))
#     colnames(wide) <- c("year", "SSB", "Total Biomass")
#     all_biom <- melt(wide, id.vars = "year")
#     all_biom2 <- cbind(all_biom,
#                        error = c(run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")], 
#                                  run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]),
#                        model = rep(name))
#     return(all_biom2)
#   }
#   
#   test_biom <- rbind(ceattle_biomass(ms_priorM1$model, "ROMS"),
#                      ceattle_biomass(run_low, "ROMS - 1"),
#                      ceattle_biomass(run_high, "ROMS + 1"))
#   
#   # Put recruitment together
#   ceattle_R <- function(run, name) {
#     R <- c(run$quantities$R)
#     error <- c(run$sdrep$sd[which(names(run$sdrep$value) == "R")])
#     R_all <- as.data.frame(cbind(year = years, 
#                                  variable = rep("Recruitment"),
#                                  value = R, 
#                                  error = error, 
#                                  model = rep(name)))
#     return(R_all)
#   }
#   R_test <- rbind(ceattle_R(ms_priorM1$model, "ROMS"),
#                   ceattle_R(run_low, "ROMS - 1"),
#                   ceattle_R(run_high, "ROMS + 1"))
#   
#   # Combine biomass & recruitment and plot
#   all_popdy <- rbind(test_biom, R_test)
#   all_popdy$year <- as.numeric(all_popdy$year)
#   all_popdy$value <- as.numeric(all_popdy$value) / 1000000  # to mt/millions
#   all_popdy$error <- as.numeric(all_popdy$error) / 1000000  # to mt/millions
#   all_popdy$variable <- factor(all_popdy$variable, labels = c("SSB (mt)", "Total Biomass (mt)", "Recruitment (millions)"))
#   
#   popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
#     geom_line(aes(linetype = model)) +
#     scale_linetype_manual(values=c("dashed", "solid", "solid"), name = "model") +
#     geom_ribbon(aes(ymin=(value-(2*error)), ymax=(value+(2*error))), alpha = 0.2, color = NA) + 
#     scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
#     scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
#     ylab(" ") +
#     labs(color = "model") +
#     facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
#     theme(strip.background = element_blank(), strip.placement = "outside")
#   
#   return(popdy_plot)
# }
# 
# test_popdy <- test_plot_popdy()
# test_popdy