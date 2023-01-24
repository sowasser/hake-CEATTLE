# Sensitivity to bioenergetics parameters considered for hake
# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230111.xlsx")

# Run CEATTLE with the values as they are in the data file
intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                                 inits = NULL, # Initial parameters = 0
                                 file = NULL, # Don't save
                                 # debug = 1, # 1 = estimate, 0 = don't estimate
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # Multispecies mode
                                 phase = "default")


### Run CEATTLE with differing bionergetics parameterizations -----------------
# Adapt weight proportions to replace those in the excel file & run CEATTLE
run_ceattle <- function(CA, Tco, Tcm, init) {
  hake_intrasp$CA <- CA
  hake_intrasp$Tco <- Tco
  hake_intrasp$Tcm <- Tcm
  ceattle <- Rceattle::fit_mod(data_list = hake_intrasp,
                               inits = init, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # Multi-species mode
                               phase = "default")
  return(ceattle)
}

# Run low-cannibalism models
run_CA2 <- run_ceattle(Tco = 8, Tcm = 10.5, CA = 0.0835, init = NULL)
run_survey <- run_ceattle(Tco = 8, Tcm = 14.5, CA = 0.167, init = NULL)

fit_CEATTLE <- function(run) {
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  jnll_summary <- as.data.frame(run$quantities$jnll_comp)
  jnll_summary$sum <- rowSums(run$quantities$jnll_comp)
  return(list(fit, jnll_summary))
}

Rceattle::plot_biomass(run_CA2, add_ci = TRUE)
Rceattle::plot_biomass(run_survey, add_ci = TRUE)


### Plot population dynamics --------------------------------------------------
bio_plot <- function() {
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
  
  test_biom <- rbind(ceattle_biomass(intrasp_run, "base model"),
                     ceattle_biomass(run_CA2, "CA / 2"),
                     ceattle_biomass(run_survey, "survey temperatures"))
  
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
  R_test <- rbind(ceattle_R(intrasp_run, "base model"),
                  ceattle_R(run_CA2, "CA / 2"),
                  ceattle_R(run_survey, "survey temperatures"))
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(test_biom, R_test)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value) / 1000000  # to mt/millions
  all_popdy$error <- as.numeric(all_popdy$error) / 1000000  # to mt/millions
  all_popdy$variable <- factor(all_popdy$variable, labels = c("SSB (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line(aes(linetype = model)) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed"), name = "model") +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    ylim(0, NA) +
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")
  
  return(popdy_plot)
}

bio_plot()
