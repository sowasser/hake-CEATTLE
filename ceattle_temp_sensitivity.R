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

hake_intrasp <- Rceattle::read_data( file = "data/hake_intrasp_221216.xlsx")

# # Run CEATTLE with the values as they are in the data file
# intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                  inits = NULL, # Initial parameters = 0
#                                  file = NULL, # Don't save
#                                  # debug = 1, # 1 = estimate, 0 = don't estimate
#                                  random_rec = FALSE, # No random recruitment
#                                  msmMode = 1, # Multispecies mode
#                                  phase = "default")

# Pull out temperature from base run & create temps for sensitivities
temp_low <- hake_intrasp$env_data$BTempC - 2
temp_high <- hake_intrasp$env_data$BTempC + 2

run_ceattle <- function(temp, df) {
  df$env_data$BTempC <- temp
  ceattle <- Rceattle::fit_mod(data_list = df,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # Multi-species mode
                               phase = "default")
  return(ceattle)
}

# Run low-cannibalism models
run_intrasp <- run_ceattle(hake_intrasp$env_data$BTempC, hake_intrasp)
run_low <- run_ceattle(temp_low, hake_intrasp)
run_high <- run_ceattle(temp_high, hake_intrasp)


# Plot population dynamics across model runs ----------------------------------
years <- 1988:2019

test_plot_popdy <- function() {
  # Pull out SSB & overall biomass from CEATTLE runs
  ceattle_biomass <- function(run, name) {
    ssb <- (c(run$quantities$biomassSSB) * 2)
    biom <- c(run$quantities$biomass)
    wide <- as.data.frame(cbind(years, ssb, biom))
    colnames(wide) <- c("year", "SSB", "Total Biomass")
    all_biom <- melt(wide, id.vars = "year")
    all_biom2 <- cbind(all_biom,
                       error = c(run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")], 
                                 run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]),
                       model = rep(name))
    return(all_biom2)
  }
  
  test_biom <- rbind(ceattle_biomass(run_intrasp, "ROMS"),
                     ceattle_biomass(run_low, "ROMS - 1"),
                     ceattle_biomass(run_high, "ROMS + 1"))
  
  # Put recruitment together
  ceattle_R <- function(run, name) {
    R <- c(run$quantities$R)
    error <- c(run$sdrep$sd[which(names(run$sdrep$value) == "R")])
    R_all <- as.data.frame(cbind(year = years, 
                                 variable = rep("Recruitment"),
                                 value = R, 
                                 error = error, 
                                 model = rep(name)))
    return(R_all)
  }
  R_test <- rbind(ceattle_R(run_intrasp, "ROMS"),
                  ceattle_R(run_low, "ROMS - 1"),
                  ceattle_R(run_high, "ROMS + 1"))
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(test_biom, R_test)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value) / 1000000  # to mt/millions
  all_popdy$error <- as.numeric(all_popdy$error) / 1000000  # to mt/millions
  all_popdy$variable <- factor(all_popdy$variable, labels = c("SSB (mt)", "Total Biomass (mt)", "Recruitment (millions)"))
  
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line(aes(linetype = model)) +
    scale_linetype_manual(values=c("dashed", "solid", "solid"), name = "model") +
    geom_ribbon(aes(ymin=(value-(2*error)), ymax=(value+(2*error))), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")
  
  return(popdy_plot)
}

test_popdy <- test_plot_popdy()
test_popdy