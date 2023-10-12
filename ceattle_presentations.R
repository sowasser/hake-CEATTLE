#' Plot results of CEATTLE runs (from run_ceattle.R), comparing the single-
#' species model and the multispecies model (cannibalism), which uses 
#' proportion of cannibalism calculated from diet database going back to 1988. 
#' Different specifications of M1 (fixed, estimated, or with a prior) are also 
#' explored.

# devtools::install_github("grantdadams/Rceattle", ref = "dev")
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent_dark.R")
theme_set(theme_sleek_transparent())

### Load in models - Rdata objects --------------------------------------------
load("models/ss_priorM1.Rdata")
load("models/ms_priorM1.Rdata")

### Plot multi-species vs. single-species vs. assessment ----------------------
start_yr <- ms_priorM1$model$data_list$styr
end_yr <- 2022
years <- start_yr:end_yr
all_yrs <- ms_priorM1$model$data_list$styr:ms_priorM1$model$data_list$projyr
hind_end <- 2019
assess_yr = "2020"
max_age <- ms_priorM1$model$data_list$nages

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

plot_models <- function(ms_run, ss_run, save_data = FALSE) {
  # Plot biomass & recruitment in comparison to no diet & assessment ----------
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
  
  biomass <- ceattle_biomass(ms_run, "CEATTLE - cannibalism")
  nodiet_biomass <- ceattle_biomass(ss_run, "CEATTLE - single-species")
  recruitment <- c(ms_run$quantities$R[, 1:length(start_yr:end_yr)])
  
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
  nodiet_R <- c(ss_run$quantities$R[, 1:length(start_yr:end_yr)])
  
  # Get recruitment from SS3 files, offset by 1 year
  ss3_R <- read.table(paste0("data/assessment/", assess_yr, "/recruitment.txt"))[-(1:((start_yr-1)-start)), ]
  ss3_R <- ss3_R[-nrow(ss3_R), ]
  
  # Put biomass together
  nodiet_biom <- ceattle_biomass(ss_run, "CEATTLE - single-species")
  
  # Combine all biomass sources together
  biom_all <- rbind(biomass, nodiet_biom, ss3_biom)
  biom_all$value <- biom_all$value / 1000000  # to Mt
  biom_all$error <- biom_all$error / 1000000  # to Mt
  
  # Put recruitment together
  R_wide <- data.frame(year = years, recruitment, nodiet_R)
  colnames(R_wide)[2:3] <- c("CEATTLE - cannibalism", "CEATTLE - single-species")
  R <- melt(R_wide, id.vars = "year")
  # Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS3 is 0)
  ss3_R <- as.data.frame(cbind(year = years, 
                               variable = rep("Assessment"), 
                               value = ss3_R[, 2],
                               error = ss3_R[, 3]))
  R_all <- rbind(cbind(R, error = c(ms_run$sdrep$sd[which(names(ms_run$sdrep$value) == "R")][1:length(start_yr:end_yr)], 
                                    ss_run$sdrep$sd[which(names(ss_run$sdrep$value) == "R")][1:length(start_yr:end_yr)])), 
                 ss3_R)
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
                           labels = c("Spawning Biomass (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  
  # Plot popdy ----------------------------------------------------------------
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line() +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.3, color = NA) + 
    scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2) +  
    scale_fill_viridis(discrete = TRUE, option = "plasma", begin = 0.2) + 
    geom_vline(xintercept = hind_end, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    ylab(" ") + xlab(" ") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 2, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside") 
  
  # Plot numbers-at-age -------------------------------------------------------
  nbyage <- extract_byage(ms_run$quantities$NByage, "CEATTLE - cannibalism", "numbers")
  nbyage_ss3_all <- read.csv(paste0("data/assessment/", assess_yr, "/nbyage.csv")) %>%
    filter(Yr >= start_yr & Yr <= end_yr)
  colnames(nbyage_ss3_all) <- c("year", "timing", c(0:max_age))
  nbyage_ss3_all <- nbyage_ss3_all[, 1:(max_age + 3)]
  
  nbyage_ss3_wide <- nbyage_ss3_all %>%
    group_by(year) %>%
    summarize_at(vars("0":as.character(max_age)), mean)
  
  nbyage_ss3 <- melt(nbyage_ss3_wide[, -2], id.vars = "year")
  nbyage_ss3 <- cbind(nbyage_ss3, 
                      rep("Stock Assessment", length(nbyage_ss3$year)))
  colnames(nbyage_ss3)[2:4] <- c("age", "numbers", "model")
  
  # Combine with nbyage from intrasp run
  nbyage <- nbyage[nbyage$year %in% 1980:2022,] ## remove extra projection years
  nbyage_all <- rbind(nbyage, nbyage_ss3)
  
  # Find mean numbers-at-age for each year
  nbyage_all$age <- as.numeric(nbyage_all$age)
  nbyage_mean <- nbyage_all %>%
    group_by(year, model) %>%
    summarize(mean = weighted.mean(age, numbers))
  nbyage_mean$year <- as.numeric(as.character(nbyage_mean$year))
  
  mean_nbyage_plot <- ggplot(nbyage_mean, aes(x = year, y = mean)) +
    geom_line(aes(color = model)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) 
  
  # Reduce down to millions for plotting
  nbyage_all$numbers <- nbyage_all$numbers / 1000000
  
  ## Difference between both models
  nbyage_diff <- cbind.data.frame(year = nbyage$year,
                                  age = as.numeric(nbyage$age),
                                  numbers = (nbyage$numbers - nbyage_ss3$numbers) / 1000000)
  nbyage_diff$model <- "CEATTLE - assessment"
  nbyage_diff$year <- factor(nbyage_diff$year)
  
  # Find midpoint between 0-1 relative to anomaly data for inform color palette
  mid <- (0 - min(nbyage_diff$numbers)) / (max(nbyage_diff$numbers) - min(nbyage_diff$numbers))
  nbyage_anomaly <- ggplot(nbyage_diff, aes(x=year, y=age)) +
    geom_point(aes(size = numbers, color = numbers)) +
    scale_color_gradientn(colours = c("#7024ac", "black", "#f3fc24"),
                          values = c(1.0, mid, 0)) +
    scale_x_discrete(breaks = c(1980, 1990, 2000, 2010, 2020)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    xlab("Year") + ylab("Age") +
    labs(size="Millions (n)", color="Millions (n)") +
    guides(size = guide_legend(override.aes = list(color="gray")))
  
  return(list(popdy = popdy_plot, nbyage_anomaly = nbyage_anomaly))
}

plots <- plot_models(ms_priorM1$model, ss_priorM1$model)
plots$popdy
plots$nbyage_anomaly

### Compare and plot natural mortality (M1 + M2) ------------------------------
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
      scale_fill_viridis(name = "M1 + M2", limits = c(0, 1.6), 
                         breaks = c(0.21, 1, 1.5),
                         option = "plasma") +
      geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
      coord_equal() +
      ylab("Age") + xlab("Year") +
      theme(panel.border = element_rect(colour = NA, fill = NA))
    mortality_plot
    
    # Min, max, mean natural mortality by age, for ages 1-5
    M_byage <- total_mortality %>%
      filter(age < 6) %>%
      group_by(age) %>%
      summarize(min = min(M1_M2), max = max(M1_M2), mean = mean(M1_M2))
    
    return(list(mortality_plot, M_byage, M1, total_mortality))
  }
}

# Cannibalism with prior on M1
ms_prior_mort <- mortality(ms_priorM1$model, type = "multi-species")
ms_prior_mort[[1]]
ms_prior_mort[[2]]

### Relative SSB / depletion --------------------------------------------------
relative_SSB <- function(model, label) {
  df <- data.frame(t(model$quantities$biomassSSB / model$quantities$SB0))
  df$year <- rownames(df)
  rownames(df) <- NULL
  df$model <- label
  df$year <- as.numeric(df$year)
  quantile <- quantile(df$Hake, probs = c(0.025, 0.5, 0.975))
  return(list(df = df, quantile = quantile))
}

relativeSSB_ss <- relative_SSB(ss_priorM1$model, "single-species")
relativeSSB_ss$quantile

relativeSSB_ms <- relative_SSB(ms_priorM1$model, "cannibalism")
relativeSSB_ms$quantile

relativeSSB_plot <- rbind(relativeSSB_ms$df, relativeSSB_ss$df) %>%
  ggplot(.) +
  geom_line(aes(x = year, y = Hake, color = factor(model))) +
  geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  scale_color_viridis(discrete = TRUE, option = "plasma", 
                      begin = 0.6, direction = -1, name = "model", 
                      guide = guide_legend(reverse = TRUE)) +
  ylab("Relative SB") +
  ylim(0, NA) +
  xlim(1980, 2022) +
  geom_hline(yintercept = 1, color = "gray40") +
  geom_hline(yintercept = 0.4, color = "gray40") +
  annotate("text", x = 1985, y = 0.4, label = "Management target", 
           vjust = -0.5, size = 2.5, color = "gray40") +
  geom_hline(yintercept = 0.1, color = "gray40") +
  annotate("text", x = 1987.5, y = 0.1, label = "Minimum stock size threshold", 
           vjust = -0.5, size = 2.5, color = "gray40") 
relativeSSB_plot

### M1 profile plots ----------------------------------------------------------
runs_ms <- list.files(path = "models/profile/ms")  # List of all model runs

ssb_all <- data.frame()
for(i in 1:length(runs_ms)) {
  load(paste0("models/profile/ms/", runs_ms[i]))
  ssb <- data.frame(t(run$quantities$biomassSSB))
  ssb$error <- (run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")] * 2)
  ssb$year <-rownames(ssb)
  rownames(ssb) <- NULL
  ssb$M1 <- round(run$quantities$M1[1, 1, 1], digits = 2)
  ssb_all <- rbind(ssb_all, ssb)
}
colnames(ssb_all)[1] <- "SSB"

ssb_all$SSB <- ssb_all$SSB / 1000000  # to Mt
ssb_all$error <- ssb_all$error / 1000000  # to Mt
ssb_all$year <- as.numeric(ssb_all$year)
ssb_all$M1 <- factor(as.character(ssb_all$M1))

ssb_all_est <- ssb_all %>% filter(M1 == 0.23)

ssb_profile_plot <- ggplot() +
  geom_line(data = ssb_all, aes(x = year, y = SSB, color = M1)) +
  geom_ribbon(data = ssb_all, 
              aes(x = year, y = SSB, ymin=(SSB - error), ymax=(SSB + error), fill = M1), 
              alpha = 0.2, color = NA) + 
  geom_line(data = ssb_all_est, aes(x = year, y = SSB), linewidth = 1, color = "white") +
  xlim(1980, 2022) +
  ylab("Spawning Biomass (Mt)") + xlab("Year") + labs(color = "M1") +
  scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2, 
                      guide = guide_legend(reverse = TRUE)) +  
  scale_fill_viridis(discrete = TRUE, option = "plasma", begin = 0.2,
                     guide = guide_legend(reverse = TRUE))  
ssb_profile_plot

# Plot fit to survey index 
init_surv <- run$data_list$srv_biom %>% filter(Year > 1)  # input survey biomass

survey_biom <- function(run, name) {
  srv <- data.frame(year = 1995:2019,
                    biomass = run$quantities$srv_bio_hat,
                    log_sd = run$quantities$srv_log_sd_hat,
                    M1 = rep(name, length(1995:2019)))
  return(srv)
}

srv_all <- data.frame()
for(i in 1:length(runs_ms)) {
  load(paste0("models/profile/ms/", runs_ms[i]))
  srv_out <- survey_biom(run = run, 
                         name = round(run$quantities$M1[1, 1, 1], digits = 2))
  srv_all <- rbind(srv_all, srv_out)
}

assess_srv <- read.csv(paste0("data/assessment/2020/survey_out.csv"))

srv_all$M1 <- factor(as.character(srv_all$M1))

survey_profile_plot <- ggplot() +
  geom_pointrange(data = init_surv, 
                  aes(x = Year, y = Observation,
                      ymin = exp(log(Observation) - 1.96*Log_sd),
                      ymax = exp(log(Observation) + 1.96*Log_sd)),
                  fatten = 5, color = "white") +
  geom_line(data = srv_all, aes(x = year, y = biomass, color = M1), alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2, 
                      guide = guide_legend(reverse = TRUE)) +  
  geom_line(data = assess_srv, aes(x = Yr, y = Exp), linetype = "dashed", color = "white") +
  xlab("year") + ylab("Index of Abundance") +
  ylim(0, NA) 
survey_profile_plot

### CEATTLE varying diet ------------------------------------------------------
# Plot diet proportion
wts <- ms_priorM1$model$data_list$UobsWtAge %>% 
  group_by(Pred_age, Prey_age) %>%
  summarize(wt_prop = mean(Stomach_proportion_by_weight))

wt05 <- scales::rescale_max(wts$wt_prop, to = c(0, 0.005))
wt10 <- scales::rescale_max(wts$wt_prop, to = c(0, 0.1))
wt50 <- scales::rescale_max(wts$wt_prop, to = c(0, 0.5))
wt75 <- scales::rescale_max(wts$wt_prop, to = c(0, 0.75))

# Look at diet proportions
prop <- as.data.frame(cbind(wts, wt05 = wt05, wt10 = wt10, wt50 = wt50, wt75 = wt75))
colnames(prop)[3:7] <- c("Base Cannibalism Model", "0.5% Cannibalism", "10% Cannibalism", 
                         "50% Cannibalism", "75% Cannibalism")
prop_all <- reshape2::melt(prop, id.vars = c("Pred_age", "Prey_age"))
prop_all$Pred_age[prop_all$Pred_age >= 15] <- "15+"
prop_all$Pred_age <- factor(as.character(prop_all$Pred_age),
                            levels = c(as.character(1:14), "15+"))
prop_all$variable <- factor(prop_all$variable, 
                            levels = c("0.5% Cannibalism", "10% Cannibalism", 
                                       "50% Cannibalism", "75% Cannibalism",
                                       "Base Cannibalism Model"))
stomach_props <- ggplot(prop_all, aes(x=Prey_age, y=value, fill=variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "plasma", begin = 0.2) + 
  ylab("stomach proportion") + xlab("prey age") +
  facet_wrap(~Pred_age, ncol = 5)
stomach_props

# Plot popdy
load("models/sensitivity/diet/run_wt05_prior.Rdata")
load("models/sensitivity/diet/run_wt10_prior.Rdata")
load("models/sensitivity/diet/run_wt50_prior.Rdata")
load("models/sensitivity/diet/run_wt75_prior.Rdata")

years <- 1980:2100
max_age <- max_age <- ms_priorM1$model$data_list$nages
models <- list(ms_priorM1$model, run_wt05_prior$model, run_wt10_prior$model, 
               run_wt50_prior$model, run_wt75_prior$model)
names <- c("Base Cannibalism Model", "0.5% Cannibalism", "10% Cannibalism", 
           "50% Cannibalism", "75% Cannibalism")

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

test_biom <- rbind(ceattle_biomass(models[[1]], names[1]),
                   ceattle_biomass(models[[2]], names[2]),
                   ceattle_biomass(models[[3]], names[3]),
                   ceattle_biomass(models[[4]], names[4]),
                   ceattle_biomass(models[[5]], names[5]))

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

R_test <- rbind(ceattle_R(models[[1]], names[1]),
                ceattle_R(models[[2]], names[2]),
                ceattle_R(models[[3]], names[3]),
                ceattle_R(models[[4]], names[4]),
                ceattle_R(models[[5]], names[5]))

plot_popdy <- function(biom, R) {
  all_popdy <- rbind(biom, R)    # Combine biomass, recruitment
  all_popdy$value <- as.numeric(all_popdy$value) / 1000000  # to mt/millions
  all_popdy$error <- as.numeric(all_popdy$error) / 1000000  # to mt/millions
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$variable <- factor(all_popdy$variable, 
                               labels = c("Spawning Biomass (Mt)", "Total Biomass (Mt)", 
                                          "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    geom_line(aes(linetype = model)) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed")) +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2) +  
    scale_fill_viridis(discrete = TRUE, option = "plasma", begin = 0.2) + 
    ylim(0, NA) + 
    xlim(1980, 2022) +
    ylab(" ") + xlab("Year") +
    labs(color = "Model", fill = "Model", linetype = "Model") +
    facet_wrap(~variable, ncol = 2, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")
  return(popdy_plot)
}
diet_sensitivity <- plot_popdy(test_biom, R_test)
diet_sensitivity

### Varying time periods ------------------------------------------------------
load("models/sensitivity/time-varying/run_90s_prior.Rdata")
load("models/sensitivity/time-varying/run_recent_prior.Rdata")

# Year ranges
all_years <- 1980:2100
no_proj <- 1980:2019

# Plot population dynamics 
timing_plot_popdy <- function(run_high, run_low, ms_model) {
  # Pull out SSB & overall biomass from CEATTLE runs
  ceattle_biomass <- function(run, name, years) {
    ssb <- (c(run$quantities$biomassSSB) * 2)
    biom <- c(run$quantities$biomass)
    wide <- as.data.frame(cbind(all_years, ssb, biom))
    colnames(wide) <- c("year", "SSB", "Total Biomass")
    all_biom <- melt(wide, id.vars = "year")
    all_biom2 <- cbind(all_biom,
                       error = c(run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")], 
                                 run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]),
                       model = rep(name)) %>%
      filter(year %in% years)
    return(all_biom2)
  }
  
  test_biom <- rbind(ceattle_biomass(ms_model, "Base Cannibalism Model", all_years),
                     ceattle_biomass(run_high, "High (1988-1999)", 1980:1999),
                     ceattle_biomass(run_low, "Low (2005-2019)", 2005:2019))
  
  # Put recruitment together
  ceattle_R <- function(run, name, years) {
    R <- c(run$quantities$R)
    error <- c(run$sdrep$sd[which(names(run$sdrep$value) == "R")])
    R_all <- as.data.frame(cbind(year = all_years, 
                                 variable = rep("Recruitment"),
                                 value = R, 
                                 error = error, 
                                 model = rep(name))) %>%
      filter(year %in% years)
    return(R_all)
  }
  R_test <- rbind(ceattle_R(ms_model, "Base Cannibalism Model", all_years),
                  ceattle_R(run_high, "High (1988-1999)", 1980:1999),
                  ceattle_R(run_low, "Low (2005-2019)", 2005:2019))
  
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(test_biom, R_test)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value) / 1000000  # to mt/millions
  all_popdy$error <- as.numeric(all_popdy$error) / 1000000  # to mt/millions
  
  # Plot population dynamics
  all_popdy$variable <- factor(all_popdy$variable, 
                               labels = c("Spawning Biomass (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  all_popdy$model <- factor(all_popdy$model, levels = c("Low (2005-2019)", "High (1988-1999)", "Base Cannibalism Model"))
  
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line() +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2) +  
    scale_fill_viridis(discrete = TRUE, option = "plasma", begin = 0.2) + 
    geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    xlim(1980, 2022) +
    ylab(" ") + xlab("Year") +
    labs(color = "Model", fill = "Model") +
    facet_wrap(~variable, ncol = 2, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")
  
  return(popdy_plot)
}

time_sensitivity <- timing_plot_popdy(run_high = run_90s_prior$model, 
                                      run_low = run_recent_prior$model,
                                      ms_model = ms_priorM1$model)

### Bioenergetics (temp-dependent consumption) --------------------------------
proxy_params <- read.csv("data/bioenergetics/proxy_bioen_params.csv")

eq_1 <- proxy_params[c(2:4), ]
eq_2 <- proxy_params[c(1, 5:6), ]
eq_2$Tco <- as.numeric(eq_2$Tco)
eq_2$Tcm <- as.numeric(eq_2$Tcm)

temp_dependent <- function(Qc, Tco, Tcm, temp) {
  z <- (log(Qc) * (Tcm - Tco)) 
  y <- (log(Qc) * (Tcm - Tco + 2))
  
  x <- (((z^2) * (1 + (1 + (40/y))^0.5)^2) / 400)
  
  consumption <- c()
  for(i in temp) { 
    V <- ((Tcm - i) / (Tcm - Tco))
    rate <- ((V^x) * exp(x * (1 - V)))
    consumption <- c(consumption, rate)
  }
  
  return(consumption)
}

# Run equation with cod, adult, juvenile pollock parameters & hake estimates
temp_range <- seq(1, 20, by = 0.1)  # Hake thermal range from survey data
spp_temp_wide <- cbind(temp_dependent(eq_2[1, 5], eq_2[1, 6], eq_2[1, 7], temp_range), 
                       temp_dependent(eq_2[2, 5], eq_2[2, 6], eq_2[2, 7], temp_range), 
                       temp_dependent(2.5, 8, 14.5, temp_range),  # Hake estimates w/ temp from acoustic series
                       temp_dependent(2.5, 8, 10.5, temp_range))  # Hake estimates w/ kriged temp
colnames(spp_temp_wide) <- c("Atlantic cod - fb4", 
                             "Pollock - fb4", 
                             "hake - survey temp",
                             "hake - kriged temp")
spp_temp <- melt(as.data.frame(spp_temp_wide))
spp_temp <- cbind(spp_temp, temp = rep(temp_range, times=4))

# # Distinguish between literature values & estimated value for Hake (when that's ready)
# spp_temp <- cbind(spp_temp, ref = c(rep("a", times = (length(temp_range) * 3)), 
#                                     rep("b", times = length(temp_range))))

# Plot consumption rate                 
temp_rate <- ggplot(spp_temp, aes(x=temp, y=value)) +
  geom_line(aes(color=variable, linetype = variable), linewidth=1) +
  # Following lines for distinguishing between lit & estimated hake values
  # geom_line(aes(color=variable, linetype=ref), size=1) +
  # scale_linetype_manual(values=c("longdash", "solid"), guide="none") +  
  scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2) +  
  xlab("temperature") + ylab("specific rate") +
  labs(color = "species", linetype = "species")
temp_rate

### Allometric mass -----------------------------------------------------------
proxy_params <- read.csv("data/bioenergetics/proxy_bioen_params.csv")

eq_1 <- proxy_params[c(2:4), ]
eq_2 <- proxy_params[c(1, 5:6), ]
eq_2$Tco <- as.numeric(eq_2$Tco)
eq_2$Tcm <- as.numeric(eq_2$Tcm)

allometric_mass <- function(CA, CB, weights) {
  Cmax <- c()
  for(i in weights) {
    mass <- (CA * (i^CB))
    Cmax <- c(Cmax, mass)
  }
  return(Cmax)
}

# Run equation with cod, adult, juvenile pollock parameters & hake estimates
wt_range <- 10:400
spp_mass_wide <- cbind(allometric_mass(eq_2[1, 3], eq_2[1, 4], wt_range), 
                       allometric_mass(eq_2[2, 3], eq_2[2, 4], wt_range), 
                       allometric_mass(0.167, -0.460, wt_range),  # hake estimate from Francis (1983)
                       allometric_mass(0.0835, -0.460, wt_range))  # Francis (1983) estimate w/ CA/2  
colnames(spp_mass_wide) <- c("Atlantic cod - fb4", 
                             "Pollock (adult) - fb4", 
                             "Hake - Francis",
                             "Hake - Francis, CA/2")
spp_mass <- melt(as.data.frame(spp_mass_wide))
spp_mass <- cbind(spp_mass, weight = rep(wt_range, times=4))

# # Distinguish between literature values & estimated value for Hake (when that's ready)
# spp_mass <- cbind(spp_mass, ref = c(rep("a", times = (length(weights) * 5)), 
#                                     rep("b", times = length(weights))))

# Plot allometric mass
mass_rate <- ggplot(spp_mass, aes(x=weight, y=value)) +
  geom_line(aes(color=variable), linewidth=1) +
  # Following lines for distinguishing between lit & estimated hake values
  # geom_line(aes(color=variable, linetype=ref), size=1) +
  # scale_linetype_manual(values=c("longdash", "solid"), guide="none") +  
  scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2) +  
  ylab("Specific Rate") + xlab("Weight (g)") +
  labs(color = "species", linetype = "species")
mass_rate

### Mean temperature ----------------------------------------------------------
survey_temp <- read.csv("data/temperature/temp_100_sophia.csv")[, -c(2:4)]
ROMS <- read.csv("data/temperature/ROMS_mean.csv")
temp_kriged <- read.csv("data/temperature/temp_100_matched_sophia.csv")

missing_years <- c(1996, 1997, 1999, 1999, 2000, 2002, 2002, 2004, 2006, 2008,2010, 
                   2014, 2016, 2018, 2020)

# Find mean from survey
survey_mean <- survey_temp %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100))

survey <- cbind(survey_mean, rep("survey", length(survey_mean$mean_temp)))
colnames(survey) <- c("year", "mean_temp", "source")

ROMS <- cbind(ROMS, rep("ROMS", length(ROMS$mean_temp)))
colnames(ROMS) <- c("year", "mean_temp", "source")

# Combine together and sort by year
CEATTLE_temp <- rbind(ROMS, survey)
CEATTLE_temp <- CEATTLE_temp[order(CEATTLE_temp$year), ]

# Only include rows where hake were found
temp_hake <- temp_kriged %>%
  filter(hake_biomass > 0)

temp_kriged_mean <- temp_kriged %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100_kriged))

temp_hake_mean <- temp_hake %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100_kriged))

temp_weighted <- temp_hake %>% group_by(year) %>%
  summarise(mean_temp = weighted.mean(temp_100_kriged, hake_biomass))

# Combine into 1 dataset with labeled data sources, then plot
means <- rbind(ROMS[, 1:2], survey_mean, temp_weighted)  # subset to model years
means <- cbind(means, c(rep("ROMS", length(1:41)),
                        rep("survey", 13), 
                        rep("kriged, biomass weighted", 12)))
colnames(means)[3] <- "dataset"
means$dataset <- factor(means$dataset, levels = c("ROMS", "survey", "kriged, biomass weighted"))

mean_temp_compared <- ggplot(means, aes(x=year, y=mean_temp)) +
  geom_point(aes(color=dataset, shape = dataset), size=2) +
  geom_line(aes(color=dataset), linewidth=1, alpha = 0.3) +
  ylim(0, NA) +
  scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2) +   
  ylab("mean temperature")
mean_temp_compared

### Diet proportion -----------------------------------------------------------
# Read in full aged dataset
aged_dataset <- read.csv("data/diet/CCTD_FEAT_combined.csv")

# Check numbers-at-age for diet data
numbers <- aged_dataset %>% 
  group_by(predator_age) %>%
  summarize(n = n())

# Set accumulation age back to 15 to deal with low sample sizes
aged_dataset$predator_age[aged_dataset$predator_age > 15] <- 15

# Create overall intraspecies predation dataset 
# Find the hake proportion for each predator
aged_wt <- aged_dataset %>%
  group_by(Predator_ID) %>%
  mutate(stomach_wt = sum(prey_wt, na.rm = TRUE)) %>%
  mutate(hake_prey_prop = if_else(prey_name == "Pacific Hake", prey_wt / stomach_wt, 0)) %>%
  select(Predator_ID, year, predator_age, prey_name, prey_age, hake_prey_prop) %>%
  distinct() # Remove duplicate rows - same pred ID, multiple hake prey 

# Total number of stomachs
stomach_n <- aged_wt %>%
  group_by(predator_age) %>%
  summarize(sample_size = n())

# Calculate average as sum of proportions per pred/prey age combo / number of stomachs per predator age
hake_prop <- aged_wt %>%
  group_by(predator_age, prey_age) %>%
  summarize(sum_prop = sum(hake_prey_prop)) %>%
  filter(!is.na(prey_age)) %>%
  left_join(stomach_n) %>%
  mutate(wt_prop = sum_prop / sample_size)

# Plot diet data
diet_plot <- ggplot(hake_prop, aes(x=as.factor(predator_age), y=wt_prop, fill=as.factor(prey_age))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis(discrete = TRUE, option = "plasma", begin = 0.2) +
  scale_x_discrete(limits = factor(1:15)) +  # add in missing predator ages
  xlab("predator hake age") + ylab("diet proportion by weight") +
  labs(fill = "prey hake age")
diet_plot

# Instances of cannibalism 
# Read in & update CCTD data from SWFSC 
CCTD_pred <- read.csv("data/diet/CCTD/hake_aged_pred.csv")
CCTD_prey <- read.csv("data/diet/CCTD/hake_aged_prey.csv")

# Combine predators and prey together, generalize prey names to "other"
CCTD_all <- merge(CCTD_pred, CCTD_prey, by = "Predator_ID") %>%
  mutate(prey_name = ifelse(Prey_Com_Name == "Pacific Hake", "Pacific Hake", "other")) %>%
  select(Predator_ID, Month, Year, Latitude, Longitude, pred_ages, prey_name, 
         Prey_Weight_g, Prey_Length1, prey_ages) %>%
  filter(Year <= 2004)  # remove years covered by FEAT 

# Add column for source of data
CCTD_all <- cbind(CCTD_all, source = rep("CCTD", length(CCTD_all[, 1])))

### Read in & update FEAT data from NWFSC 
FEAT_data <- read.csv("data/diet/FEAT/all_FEAT_diet.csv")
FEAT_prey_lengths <- read.csv("data/diet/FEAT/hake_prey_lengths_FEAT.csv")

# Add prey lengths based on stomach ID
FEAT_all <- merge(FEAT_data, FEAT_prey_lengths, by = "stomach_uuid", all.x = TRUE)

# Expand to month and year columns
FEAT_all$tow_timestamp <- as.Date(FEAT_all$tow_timestamp)
FEAT_all$month <- lubridate::month(lubridate::ymd(FEAT_all$tow_timestamp))
FEAT_all$year <- lubridate::year(lubridate::ymd(FEAT_all$tow_timestamp))

# Generalize prey categories to Pacific Hake & "other", select columns
FEAT_all <- FEAT_all %>%
  mutate(prey_name = ifelse(prey_category_long == "Gadiformes", "Pacific Hake", "other")) %>%
  select(stomach_uuid, month, year, tow_latitude, tow_longitude, predator_age, 
         prey_name, content_wt_g, measure_value)

# Switch to age 20 accumulator age
FEAT_all <- FEAT_all %>% 
  filter(predator_age != "(blank)") %>%
  mutate(predator_age = as.numeric(predator_age)) %>%
  mutate(predator_age = ifelse(predator_age > 20, 20, predator_age))

#Combine datasets
# Add column for prey ages. All are age 1 according to calculations in the 
# diet_aging.R script
FEAT_all <- cbind(FEAT_all, 
                  prey_ages = rep(1, length(FEAT_all[, 1])),
                  source = rep("FEAT", length(FEAT_all[, 1])))

labels <- c("Predator_ID", "month", "year", "latitude", "longitude", "predator_age", 
            "prey_name", "prey_wt", "prey_length", "prey_age", "source")

colnames(CCTD_all) <- labels
colnames(FEAT_all) <- labels

all_data <- rbind(CCTD_all, FEAT_all)

# Remove all prey ages for prey that isn't hake & remove 1980 (year w/ no ages)
all_data$prey_age[all_data$prey_name == "other"] <- NA

hake_data <- all_data %>%
  filter(prey_name == "Pacific Hake") 

# Look at instances of no age - all weights within year 1 range except for 1
no_age <- hake_data %>% filter(is.na(prey_length))
# Replace NAs with age 1
all_data$prey_age[is.na(all_data$prey_age) & all_data$prey_name == "Pacific Hake"] <- 1 

# Add estimate of correct age for very big hake prey
all_data$prey_age[all_data$prey_wt > 300 & all_data$prey_name == "Pacific Hake"] <- 5

# Rewmove 1980 w/ no ages
all_data <- all_data %>% filter(year != 1980)

# Plot cannibalism rate
loc_n <- all_data %>%
  group_by(year, latitude, longitude) %>%
  summarize(n_all = n()) %>%
  filter(!is.na(year))

loc_n_cannibalism <- all_data %>%
  filter(prey_name == "Pacific Hake") %>%
  group_by(year, latitude, longitude) %>%
  summarize(n_cannibalism = n()) %>%
  filter(!is.na(year))

loc_n_all <- left_join(loc_n, loc_n_cannibalism)
loc_n_all$n_cannibalism[is.na(loc_n_all$n_cannibalism)] <- 0
loc_n_all$prop <- loc_n_all$n_cannibalism / loc_n_all$n_all
loc_n_all <- loc_n_all %>% arrange(prop)

# New plot of cannibalism rate by year
annual_rate <- loc_n_all %>%
  group_by(year) %>%
  summarize(n = sum(n_all),
            prop = mean(prop)) %>%
  ggplot(., aes(x = year, y = n, fill = prop)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(option = "plasma", limits = c(0, 0.4), begin = 0.2) +
  ylab("Stomachs (n)") +
  labs(fill = "Cannibalism Rate")
annual_rate

### Compare original proportions to weighted proportions ----------------------
original <- read.csv("data/diet/diet_for_CEATTLE_original.csv")

# Read in analyzed data if just fixing plots
ceattle_all <- read.csv("data/diet/Dirichlet/Dirichlet_all_years.csv")
ceattle_90s <- read.csv("data/diet/Dirichlet/Dirichlet_90s.csv")
ceattle_recent <- read.csv("data/diet/Dirichlet/Dirichlet_recent.csv")

new_df <- function(df, name) {
  new <- cbind(pred_age = df$Pred_age, 
               prey_age = df$Prey_age, 
               prop = df$Stomach_proportion_by_weight, 
               data = rep(name))
  return(new)
} 

# Neater comparison graph for write-up
comp2 <- rbind(new_df(ceattle_all, "All Years"),
               new_df(ceattle_90s, "High (1991-1999)"),
               new_df(ceattle_recent, "Low (2005-2019)"))

comp2 <- as.data.frame(comp2) 
comp2$pred_age <- as.factor(comp2$pred_age)
comp2$prey_age <- as.numeric(comp2$prey_age)
comp2$prop <- as.numeric(comp2$prop)
comp2$data <- factor(comp2$data, 
                     levels = c("All Years", "High (1991-1999)", "Low (2005-2019)"))
comp2 <- comp2 %>% filter(prey_age <= 5)

comp2_plot <- ggplot(comp2, aes(x=pred_age, y=prop, fill=factor(prey_age))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits = factor(1:15), 
                   breaks = c(1, 5, 10, 15), 
                   labels = c("1", "5", "10", "15+")) +  # add in missing predator ages
  scale_y_continuous(limits = c(0, 1), labels = scales::label_number(accuracy = NULL)) +
  scale_fill_viridis(discrete = TRUE, option = "plasma", begin = 0.2) +
  xlab("Predator Hake Age") + ylab("Diet Proportion") +
  labs(fill = "Prey Hake Age") +
  facet_wrap(~ data)
comp2_plot

### Plot timing of sample collection ----------------------------------------
time_n <- all_data %>%
  group_by(year, month) %>%
  summarize(n_all = n()) %>%
  filter(!is.na(year))

time_n_cannibalism <- all_data %>%
  filter(prey_name == "Pacific Hake") %>%
  group_by(year, month) %>%
  summarize(n_cannibalism = n()) %>%
  filter(!is.na(year)) 

time_n_all <- left_join(time_n, time_n_cannibalism)
time_n_all$n_cannibalism[is.na(time_n_all$n_cannibalism)] <- 0
time_n_all$prop <- time_n_all$n_cannibalism / time_n_all$n_all
time_n_all <- time_n_all %>% arrange(year)
# time_n_all$month <- factor(time_n_all$month)

data_years <- c(1988:1992, 1995:1999, 2002, 2005, 2007, 2009, 2011:2013, 2015, 2017, 2019)
all_months <- cbind.data.frame(year = rep(data_years, each = 12),
                               month = rep(1:12, times = length(data_years)))
time_all_months <- left_join(all_months, time_n_all)
time_all_months[is.na(time_all_months)] <- 0

overall_monthly <- time_all_months %>% 
  group_by(month) %>%
  summarize(n_all = sum(n_all), n_cannibalism = sum(n_cannibalism)) 
overall_monthly$prop <- overall_monthly$n_cannibalism / overall_monthly$n_all
overall_monthly$prop[is.nan(overall_monthly$prop)] <- 0

monthly_rate <- ggplot(overall_monthly, aes(x = month, y = n_all, fill = prop)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = factor(1:12), breaks = c(2, 4, 6, 8, 10, 12)) +
  scale_fill_viridis(option = "plasma", limits = c(0, 0.52), begin = 0.2) +
  xlab("Sampling Month") + ylab("Stomachs (n)") +
  labs(fill = "Cannibalism Rate") 
monthly_rate

### Location and timing plots -------------------------------------------------
# Create a plot of location of observations by latitude and longitude
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry

locations <- ggplot(data = world) +
  geom_sf(fill = "gray30") +
  geom_point(data = loc_n_all, aes(x = longitude, y = latitude, color = prop, size = n_all)) +
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) +
  scale_x_continuous(breaks = seq(-135, -115, by = 5)) +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) +
  scale_color_viridis(option = "plasma", begin = 0.2) +
  xlab(" ") + ylab(" ") + labs(color = "Cannibalism Rate", size = "Stomachs (n)") +
  guides(size = guide_legend(override.aes = list(color="gray"))) +
  facet_wrap(~year, ncol = 7)

# Inset timing plots in yearly location plots 
# Tutorial here: https://www.blopig.com/blog/2019/08/combining-inset-plots-with-facets-using-ggplot2/
get_inset <- function(df) {
  # Create plot for the inset 
  plot <- ggplot(df, aes(x = factor(month), y = n_all, fill = prop)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(limits = factor(1:12), breaks = c(1, 12), labels = c("Jan", "Dec")) +
    scale_fill_viridis(option = "plasma", limits = c(0, 1), begin = 0.2) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=rel(0.9)),  # inset axis tick font size
          plot.background = element_rect(fill='transparent', color=NA)) + # transparent so no overlap w/map
    theme(legend.position="none") 
  return(plot)
}

# Function for defining how the inset will be positioned
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

inset_plot <- get_inset(time_n_all)  # Actually create insets

# How the insets will be mapped on to the main plots (applying above function)
insets <- time_all_months %>% 
  split(f = .$year) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.)), 
    data = data.frame(year=unique(.$year)),
    ymin = 30.7, ymax = 40, xmin = -141, xmax = -124))  # position of insets

# Bring everything together - add insets on to main plot (locations, created above)
location_timing <- locations +
  coord_sf(xlim = c(-140, -115), ylim = c(31, 56), expand = FALSE) + 
  scale_x_continuous(breaks = c(-135, -120)) +
  scale_y_continuous(breaks = c(32, 55)) +
  insets

### Hake predators ------------------------------------------------------------
path <- "data/diet/CCTD/v4/"

# Read in all predator and prey data
predator <- read.csv(paste0(path, "predator_information_v4.csv"))
prey_comp <- read.csv(paste0(path, "prey_composition_v4.csv"))

pred_prey <- merge(predator, prey_comp, by = "Predator_ID")

# Calculate highest predation by relative occurrence
high_n <- pred_prey %>%
  group_by(Predator_Com_Name, Prey_Com_Name) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(Prey_Com_Name == "Pacific Hake") %>%
  arrange(-freq)

# Plot highest predators
hake_pred_plot <- ggplot(high_n, aes(x = reorder(Predator_Com_Name, freq), y = freq)) +
  geom_bar(position = "dodge", stat = "identity", show.legend = FALSE, fill = "#6921a8") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  xlab(" ") + ylab("Relative frequency of hake predation") 
hake_pred_plot

### Save plots (when not experimenting) ---------------------------------------
ggsave(filename="plots/presentations/popdyn_M1prior.png", plots$popdy, 
       width=200, height=90, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/nbyage_anomaly.png", plots$nbyage_anomaly, 
       width=200, height=90, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/M.png", ms_prior_mort[[1]], 
       width = 160, height = 70, units = "mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/relative_SSB.png", relativeSSB_plot, 
       width=150, height=80, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/M1_profile_SSB.png", ssb_profile_plot, 
       width=140, height=80, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/M1_profile_survey.png", survey_profile_plot, 
       width=140, height=80, units="mm", dpi=300, bg = "transparent")
ggsave(filename = "plots/presentations/sensitivity_prop.png", stomach_props,
       width=270, height=140, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/diet_sensitivity.png", diet_sensitivity, 
       width=200, height=90, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/time_sensitivity.png", time_sensitivity, 
       width=200, height=90, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/temp_consumption.png", temp_rate,
       bg = "transparent", width=160, height=90, units="mm", dpi=300)
ggsave(filename="plots/presentations/allometric_mass.png", mass_rate,
       bg = "transparent", width=160, height=90, units="mm", dpi=300)
ggsave(filename="plots/presentations/mean_temp_compared.png", mean_temp_compared,
       bg = "transparent", width=180, height=90, units="mm", dpi=300)
ggsave(filename = "plots/presentations/cannibalism_overall.png", diet_plot, 
       bg = "transparent", width=200, height=120, units="mm", dpi=300)
ggsave(filename = "plots/presentations/hake_cannibalism_rate.png", annual_rate,
       bg = "transparent", width=190, height=100, units="mm", dpi=300)
ggsave(filename = "plots/presentations/Dirichlet_comp_pretty.png", comp2_plot, 
       bg = "transparent", width=170, height=50, units="mm", dpi=300)
ggsave(filename = "plots/presentations/monthly_rate.png", monthly_rate,
       bg = "transparent", width=190, height=100, units="mm", dpi=300)
ggsave(filename = "plots/presentations/locations_timing.png", location_timing,
       bg = "transparent", width=300, height=200, units="mm", dpi=300)
ggsave(filename = "plots/presentations/hake_predators.png", hake_pred_plot, 
       bg = "transparent", width=170, height=100, units="mm", dpi=300)
