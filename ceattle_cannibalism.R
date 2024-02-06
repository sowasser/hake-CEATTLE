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
library(ggview)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

### Load in models - Rdata objects --------------------------------------------
# load("models/ss_fixM1.Rdata")
# load("models/ss_estM1.Rdata")
load("models/ss_priorM1.Rdata")
# load("models/ms_estM1.Rdata")
# load("models/ms_fixM1.Rdata")
load("models/ms_priorM1.Rdata")

# Compare fits
model_fits <- rbind(# cbind(model = "SS est M1", ss_estM1$fit),
                    # cbind(model = "SS fix M1", ss_fixM1$fit),
                    cbind(model = "SS prior M1", ss_priorM1$fit),
                    # cbind(model = "MS est M1", ms_estM1$fit),
                    # cbind(model = "MS fix M1", ms_fixM1$fit),
                    cbind(model = "MS prior M1", ms_priorM1$fit))

model_summary <- list(ss_priorM1$summary, ms_priorM1$summary) %>%
  reduce(full_join, by = "component") 
model_summary$NLL.x <- round(model_summary$NLL.x, digits = 1)
model_summary$NLL.y <- round(model_summary$NLL.y, digits = 1)
colnames(model_summary) <- c("component", "SS prior M1", "MS prior M1")


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
  # ms_run <- ms_priorM1$model
  # ss_run <- ss_priorM1$model
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
                                  model = "2020 Assessment"))
  colnames(ss3_biom)[2:3] <- c("value", "error")
  ss3_biom <- ss3_biom[, c(1, 4, 2, 3, 5)]
  
  # Pull out recruitment
  nodiet_R <- c(ss_run$quantities$R[, 1:length(start_yr:end_yr)])
  
  # Get recruitment from SS3 files, offset by 1 year
  ss3_R <- read.table(paste0("data/assessment/", assess_yr, "/recruitment.txt"))[-(1:((start_yr-1)-start)), ]
  ss3_R <- ss3_R[-nrow(ss3_R), ]
  
  # Find mean difference between the model runs (for hindcast) ----------------
  rel_change <- function(df1, df2, title) {
    mean_out <- round(mean((df1 / 1000000) - (df2 / 1000000)), 3)
    SEM <- round(sd((df1 / 1000000) - (df2 / 1000000)) / sqrt(length(range)), 3)
    percent <- round(mean(((df1 - df2) / df1) * 100), 3)
    return(c(title, mean_out, SEM, percent))
  }
  
  rechange_all <- rbind(rel_change(biomass[biomass$type == "SSB" & biomass$year <= hind_end,]$value,
                                   nodiet_biomass[nodiet_biomass$type == "SSB" & biomass$year <= hind_end,]$value,
                                   "cannibalism - no diet, SSB"),
                        rel_change(biomass[biomass$type == "Total Biomass" & biomass$year <= hind_end,]$value,
                                   nodiet_biomass[nodiet_biomass$type == "Total Biomass" & biomass$year <= hind_end,]$value,
                                   "cannibalism - no diet, Total"),
                        rel_change(recruitment[1:length(start_yr:hind_end)], 
                                   nodiet_R[1:length(start_yr:hind_end)],
                                   "cannibalism - no diet, R"),
                        rel_change(nodiet_biomass[nodiet_biomass$type == "SSB" & nodiet_biomass$year <= hind_end,]$value,
                                   ss3_biom[ss3_biom$type == "SSB" & ss3_biom$year <= hind_end,]$value,
                                   "no diet - SS3, SSB"),
                        rel_change(nodiet_biomass[nodiet_biomass$type == "Total Biomass" & biomass$year <= hind_end,]$value,
                                   ss3_biom[ss3_biom$type == "Total Biomass" & ss3_biom$year <= hind_end,]$value,
                                   "no diet - SS3, Total"),
                        rel_change(nodiet_R[1:length(start_yr:hind_end)],
                                   ss3_R[1:length(start_yr:hind_end), 2],
                                   "no diet - SS3, R"))
  colnames(rechange_all) <- c("model", "mean difference", "SEM", "percent")
  
  # Plot biomass & recruitment ------------------------------------------------
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
                               variable = "2020 Assessment", 
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
  
  # Relative SSB / depletion ----------------------------------------------------
  dynB0_ss3 <- as.data.frame(read.table(paste0("data/assessment/", assess_yr, "/dyn_B0.txt")))[, 2]
  dynB0_ss3 <- mean(dynB0_ss3) / 1000000
  relSSB <- rbind.data.frame(cbind.data.frame(depletion = t(ss_run$quantities$depletionSSB)[1:43],
                                              year = years, 
                                              model = "CEATTLE - single-species"),
                             cbind.data.frame(depletion = t(ms_run$quantities$depletionSSB)[1:43],
                                              year = years, 
                                              model = "CEATTLE - cannibalism"),
                             cbind.data.frame(depletion = all_popdy[all_popdy$type == "SSB" & 
                                                                      all_popdy$model == "2020 Assessment", ]$value / (dynB0_ss3),
                                              year = years,
                                              model = "2020 Assessment"))
    
  # Set columns to match popdy dataframe
  relSSB$error <- 0
  relSSB$type <- "Relative SB"
  relSSB <- relSSB[, c(2, 5, 1, 4, 3)]
  colnames(relSSB)[3] <- "value"
  
  # Combine and set up for plotting
  all_popdy <- rbind.data.frame(all_popdy, relSSB)
  all_popdy$model <- factor(all_popdy$model, 
                            levels = c("2020 Assessment", 
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
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    ylab(" ") + xlab("Year") +
    labs(color = "Model", fill = "Model") +
    facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside") 
  
  # Plot ratio of SSB:Biomass to look for skewness in age composition ---------
  ratio <- as.data.frame(cbind(year = years,
                               assessment = ss3_ssb[, 1]/ss3_biomass[, 1],
                               cannibalism = biomass[1:length(years), 3] / 
                                 biomass[(length(years)+1):(length(years)*2), 3],
                               single_species = nodiet_biom[1:length(years), 3] / 
                                 nodiet_biom[(length(years)+1):(length(years)*2), 3]))
  
  colnames(ratio) <- c("year", "2020 Assessment", "CEATTLE - cannibalism", 
                       "CEATTLE - single-species")
  ratio2 <- melt(ratio, id.vars = "year", variable.name = "model")
  ratio2$model <- factor(ratio2$model, levels = c("Assessment", 
                                                  "CEATTLE - single-species", 
                                                  "CEATTLE - cannibalism"))
  
  ratio_plot <- ggplot(ratio2, aes(x=year, y=value, color=model)) +
    geom_line() +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylab("SSB/Biomass")
  
  # Plot the difference between model runs ------------------------------------
  ceattle_intrasp <- all_popdy %>% filter(model == "CEATTLE - cannibalism")
  ceattle_nodiet <- all_popdy %>% filter(model == "CEATTLE - single-species")
  assessment <- all_popdy %>% filter(model == "2020 Assessment")
  
  diff_intrasp <- rbind(cbind.data.frame(year = years,
                                         type = ceattle_intrasp$type,
                                         difference = ceattle_intrasp$value - ceattle_nodiet$value,
                                         models = "cannibalism - single-species"),
                        cbind.data.frame(year = years,
                                         type = ceattle_nodiet$type,
                                         difference = ceattle_nodiet$value - assessment$value,
                                         models = "single-species - assessment"))
  
  diff_plot <- ggplot(diff_intrasp, aes(x=year, y=difference, color=models)) +
    geom_line() +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    geom_vline(xintercept = 2020, linetype = 2, colour = "gray")  # Add line at end of hindcast
  
  # Plot numbers-at-age -------------------------------------------------------
  nbyage <- extract_byage(ms_run$quantities$NByage, "CEATTLE - cannibalism", "numbers")
  
  # # Read in data from no diet CEATTLE run
  # nbyage_nodiet <- extract_nbyage(ss_run, "CEATTLE - single species")

  # Read in data from SS3 & average beginning & middle of the year
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
  
  # Plot yearly nbyage
  nbyage_plot <- ggplot(nbyage_all, aes(x=year, y=age)) +
    geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
    scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
    scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
    scale_x_discrete(breaks = c(1980, 1990, 2000, 2010, 2020)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    geom_hline(yintercept = 2.5, color = "gray") +  # line at maturity ogive
    xlab("Year") + ylab("Age") + 
    labs(fill="Millions", size="Millions", color="Millions") +
    facet_wrap(~model, ncol=1)
  
  # Numbers-at-age anomaly ----------------------------------------------------
  nbyage_diff <- cbind.data.frame(year = nbyage$year,
                                  age = as.numeric(nbyage$age),
                                  numbers = (nbyage$numbers - nbyage_ss3$numbers) / 1000000)
  nbyage_diff$model <- "CEATTLE - assessment"
  nbyage_diff$year <- factor(nbyage_diff$year)
  
  # Find midpoint between 0-1 relative to anomaly data for inform color palette
  mid <- (0 - min(nbyage_diff$numbers)) / (max(nbyage_diff$numbers) - min(nbyage_diff$numbers))
  nbyage_anomaly <- ggplot(nbyage_diff, aes(x=year, y=age)) +
    geom_point(aes(size = numbers, color = numbers)) +
    scale_color_gradientn(colours = c("#5f4187", "white", "#bbdf27"),
                          values = c(1.0, mid, 0)) +
    scale_x_discrete(breaks = c(1980, 1990, 2000, 2010, 2020)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    geom_hline(yintercept = 2.5, color = "gray") +  # line at maturity ogive
    xlab("Year") + ylab("Age") +
    labs(size="Number \n(Millions)", color="Number \n(Millions)")

  # Plot comparison to survey index -------------------------------------------
  init_surv <- ms_run$data_list$srv_biom %>% filter(Year > 1)  # input survey biomass
  
  survey_biom <- function(run, name) {
    srv <- data.frame(year = 1995:hind_end,
                      biomass = run$quantities$srv_bio_hat,
                      log_sd = run$quantities$srv_log_sd_hat,
                      model = rep(name, length(1995:hind_end)))
    return(srv)
  }
  
  intrasp_srv <- survey_biom(ms_run, "CEATTLE - cannibalism")
  nodiet_srv <- survey_biom(ss_run, "CEATTLE - single-species")
  
  survey <- read.csv(paste0("data/assessment/", assess_yr, "/survey_simple.csv"))
  colnames(survey) <- c("year", "biomass", "log_sd")
  survey <- cbind(survey, model = "2020 Assessment")
  
  survey_all <- rbind(intrasp_srv, nodiet_srv, survey)
  survey_all$model <- factor(survey_all$model, levels = c("2020 Assessment", "CEATTLE - single-species", "CEATTLE - cannibalism"))
  survey_all$Observation <- survey_all
  
  survey_plot <- ggplot() +
    geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    geom_pointrange(data = init_surv, 
                    aes(x = Year, y = Observation / 1000000,
                        ymin = exp(log(Observation / 1000000) - 1.96*Log_sd),
                        ymax = exp(log(Observation / 1000000) + 1.96*Log_sd)),
                    fatten = 5) +
    geom_line(data = survey_all, aes(x = year, y = (biomass / 1000000), color = model), alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    ylim(0, NA) +
    xlab("year") + ylab("Index of Abundance")
  
  # Suitability ---------------------------------------------------------------
  suitability <- as.data.frame(as.table(ms_run$quantities$suit_main)) %>%
    filter(Var1 == "A" & Var2 == "A")
  
  suitability <- suitability[, 3:6]
  suitability$Var3 <- as.integer(factor(suitability$Var3, labels = c(1:max_age)))
  suitability$Var4 <- as.integer(factor(suitability$Var4, labels = c(1:max_age)))
  suitability$Var5 <- factor(suitability$Var5, labels = all_yrs)
  
  colnames(suitability) <- c("pred_age", "prey_age", "year", "value")
  suitability <- suitability %>% filter(prey_age < 6)
  suitability$prey_age <- factor(suitability$prey_age)
  suitability$year <- as.numeric(as.character(suitability$year))
  suitability <- suitability %>% 
    filter(year <= 2022) %>%
    group_by(pred_age, prey_age) %>%
    summarize(value = mean(value)) %>%
    filter(pred_age <= 15)
  suitability$pred_age[suitability$pred_age == 15] <- "15+"
  suitability$pred_age <- factor(suitability$pred_age,
                                 levels = c(1:14, "15+"))
  
  suit_plot <- ggplot(suitability, aes(x = pred_age, y = value, fill = prey_age)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, name = "Prey Age") +
    xlab("Predator Age") + ylab("Suitability") 
  
  # Biomass-at-age for each model run -----------------------------------------
  biombyage <- extract_byage(ms_run$quantities$biomassByage, 
                             "CEATTLE - cannibalism", "biomass")
  
  biombyage$age <- as.integer(biombyage$age)
  biombyage <- biombyage[as.numeric(as.character(biombyage$year)) <= end_yr,]
  biombyage$biomass <- biombyage$biomass / 1000000
  
  # Plot yearly biomass by age
  biombyage_plot <- ggplot(biombyage, aes(x=year, y=age)) +
    geom_point(aes(size = biomass, color = biomass)) +
    scale_color_viridis(direction = -1) +
    scale_x_discrete(breaks = c(1980, 1990, 2000, 2010, 2020)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    geom_hline(yintercept = 2.5, color = "gray") +  # line at maturity ogive
    xlab("Year") + ylab("Age") +
    labs(size="Biomass \n(Mt)", color="Biomass \n(Mt)")
  
  #Plot realized consumption --------------------------------------------------
  # Extract biomass consumed as prey
  b_consumed <- extract_byage(ms_run$quantities$B_eaten_as_prey, 
                              "CEATTLE - cannibalism", "biomass")[, -4]
  b_consumed$age <- as.integer(b_consumed$age)
  
  # Filter to only ages below 6 (the ages of hake consumed)
  b_consumed <- b_consumed %>% filter(age < 6)
  b_consumed$age <- as.factor(b_consumed$age)
  b_consumed <- b_consumed[as.numeric(as.character(b_consumed$year)) <= end_yr,]
  
  # Plot yearly biomass consumed by age
  b_consumed_plot <- ggplot(b_consumed, aes(x=year, y=(biomass / 1000000), fill = age)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    scale_x_discrete(breaks = seq(start_yr, end_yr, 3)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylab("Biomass consumed (Mt)")
  b_consumed_plot
  
  # Plot ratio of biomass consumed to approximation of predator biomass (SSB)
  # Add ages together
  yearly_consumed <- b_consumed %>%
    group_by(year) %>%
    summarize(total_biomass = sum(biomass)) %>%
    ungroup()
  
  ssb <- biomass %>% filter(type == "SSB")
  
  yearly_consumed$total_biomass <- (yearly_consumed$total_biomass / (ssb$value))
  yearly_consumed$year <- as.numeric(as.character(yearly_consumed$year))
  yearly_b_plot <- ggplot(yearly_consumed, aes(x = year, y = total_biomass)) +
    geom_line() +
    geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylab("biomass of prey / SSB")
  
  if(save_data == TRUE) {
    write.csv(biomass, paste0("data/CEATTLE/", assess_yr, "/ceattle_intrasp_biomass.csv"), row.names = FALSE)
    write.csv(nodiet_biomass, paste0("data/CEATTLE/", assess_yr, "/ceattle_nodiet_biomass.csv"), row.names = FALSE)
    write.csv(recruitment, paste0("data/CEATTLE/", assess_yr, "/ceattle_intrasp_R.csv"), row.names = FALSE)
    write.csv(nbyage, paste0("data/CEATTLE/", assess_yr, "/ceattle_intrasp_nbyage.csv"), row.names = FALSE)
  }
  
  return(list(relative_change = rechange_all, popdy = popdy_plot, ratio = ratio_plot, 
              pop_diff = diff_plot, nbyage = nbyage_plot, 
              nbyage_anomaly = nbyage_anomaly, survey = survey_plot, 
              suit = suit_plot, biombyage = biombyage_plot, 
              b_consumed = b_consumed_plot, yearly_b = yearly_b_plot))
}

plots <- plot_models(ms_priorM1$model, ss_priorM1$model)
relative_change <- plots$relative_change
plots$popdy
plots$ratio
plots$pop_diff
plots$nbyage
plots$nbyage_anomaly
plots$survey
plots$suit
plots$biombyage
plots$b_consumed
plots$yearly_b

# Combine nbyage plots togehter -----------------------------------------------
biom_n_byage <- cowplot::plot_grid(plots$biombyage, plots$nbyage_anomaly, ncol = 1, labels = c("A", "B"))

# # Plot with fixed M1
# plots_M1fixed <- plot_models(ms_fixM1$model, ss_fixM1$model)
# relative_change_M1fix <- plots_M1fixed$relative_change
# plots_M1fixed$popdy
# 
# # Plot with M1 estimated
# plots_M1est <- plot_models(ms_estM1$model, ss_estM1$model)
# relative_change_M1est <- plots_M1est$relative_change
# plots_M1est$popdy

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
    
    return(list(mortality_plot, M_byage, M1, total_mortality))
  }
}

# # Cannibalism with estimated M1
# ms_est_mort <- mortality(ms_estM1$model, type = "multi-species")
# ms_est_mort[[1]]
# ms_est_mort[[2]]
# ms_est_M1 <- ms_est_mort[[3]]
# ms_est_totM <- ms_est_mort[[4]] %>% 
#   group_by(year) %>%
#   summarize(M1_M2 = sum(M1_M2))
# 
# # Cannibalism with fixed M1
# ms_fixed_mort <- mortality(ms_fixM1$model, type = "multi-species")
# ms_fixed_mort[[1]]
# ms_fixed_mort[[2]]
# ms_fix_totM <- ms_fixed_mort[[4]] %>% 
#   group_by(year) %>%
#   summarize(M1_M2 = sum(M1_M2))

# Cannibalism with prior on M1
ms_prior_mort <- mortality(ms_priorM1$model, type = "multi-species")
ms_prior_mort[[1]]
ms_prior_mort[[2]]
ms_prior_M1 <- ms_prior_mort[[3]]
ms_prior_totM <- ms_prior_mort[[4]] %>% 
  group_by(age, year) %>%
  summarize(M1_M2 = sum(M1_M2))

ms_totM_age1 <- ms_prior_totM %>%
  filter(age == 1)

# Create dataframe for a table of total mortality (rounded and only up to age 6)
m_tot <- ms_prior_totM %>% 
  filter(age <= 6) %>%
  mutate(M1_M2 = round(M1_M2, digits = 2)) %>%
  dcast(., year ~ age)

write.csv(m_tot, file = "data/CEATTLE/total_M.csv", row.names = FALSE)

# # Single-species with estimated M1
# ss_M1 <- mortality(ss_estM1$model, type = "single-species")
# ss_prior_M1 <- mortality(ss_priorM1$model, type = "single-species")

### Reference points ----------------------------------------------------------
# Fishing mortality
eff <- rbind(data.frame(year = years,
                        eff = ms_priorM1$model$quantities$F_spp[1:length(start_yr:end_yr)],
                        model = "Cannibalism"),
             data.frame(year = years,
                        eff = ss_priorM1$model$quantities$F_spp[1:length(start_yr:end_yr)],
                        model = "Single-species"))
eff$model <- factor(eff$model, levels = c("Single-species", "Cannibalism"))

eff_plot <- ggplot(eff, aes(x=year, y=eff, color=model)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.5) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.5) +
  geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  xlab("Year") + ylab("Fishing Mortality") + labs(color = "CEATTLE Model")
eff_plot

# ms_run_Fspr <- Rceattle::fit_mod(data_list = ms_priorM1$model$data_list,
#                                  inits =  ms_priorM1$model$estimated_params, # Initial parameters from ss_run
#                                  estimateMode = 2, # Run projection only
#                                  HCR = build_hcr(HCR = 4, # NEFMC HCR
#                                                  FsprTarget = 0.4, # 0.75 * F40%
#                                                  FsprLimit = 0.4, # F40%
#                                                  Fmult = 0.75,
#                                                  Plimit = 0.2),
#                                  M1Fun = Rceattle::build_M1(M1_model = 1,
#                                                             updateM1 = TRUE,
#                                                             M1_use_prior = TRUE,
#                                                             M1_prior_mean = 0.2,
#                                                             M1_prior_sd = .1),
#                                  msmMode = 1, # Single species mode
#                                  verbose = 1)
# 
# mean(ms_run_Fspr$quantities$SB0)
# ms_run_Fspr$quantities$B0 #B0
# ms_run_Fspr$quantities$SB0 #SB0
# ms_run_Fspr$quantities$Flimit #F that gives you SPR40%
# ms_run_Fspr$quantities$SPRlimit #SPR40%

brp_comparison <- function(model, model_name) {
  df <- data.frame(c((model$quantities$B0[121] / 1000000),
                     (model$quantities$SB0[121] / 1000000),
                     (model$quantities$R0 / 1000000),
                     model$quantities$Flimit,
                     model$quantities$SPRlimit,
                     model$quantities$SPRtarget))
  row.names(df) <- c("B0", "SB0", "R0", "Flimit", "SPRlimit", "SPRtarget")
  colnames(df) <- model_name
  df <- round(df, 2)
  return(df)
}

brps <- cbind(
  # brp_comparison(model = ss_estM1$model, model_name = "SS est M1"),
  # brp_comparison(model = ss_fixM1$model, model_name = "SS fix M1"),
  brp_comparison(model = ss_priorM1$model, model_name = "SS prior M1"),
  # brp_comparison(model = ms_estM1$model, model_name = "MS est M1"),
  # brp_comparison(model = ms_fixM1$model, model_name = "MS fix M1"),
  brp_comparison(model = ms_priorM1$model, model_name = "MS prior M1")
)

data.frame(t(ms_priorM1$model$quantities$biomassSSB / ms_priorM1$model$quantities$SB0))
data.frame(t(ss_priorM1$model$quantities$biomassSSB / ss_priorM1$model$quantities$SB0))

# Plot of empirical weight-at-age data from the assessment --------------------
# Extract index used for derived quantities
weight_out <- ss_priorM1$model$data_list$wt %>% filter(Wt_index == 2)  
weight_out <- weight_out[, 5:20]  # only keep needed columns & ages 1-15
# Reshape dataframe
colnames(weight_out) <- c("Year", 1:15)
weight <- melt(weight_out, id.vars = "Year", variable.name = "Age",
               value.name = "Weight") 
# Find mean and combine with dataframe @ 2020 so it can be plotted alongside
weight_mean <- weight %>% 
  group_by(Age) %>%
  summarize(mean = mean(Weight))
weight_mean <- cbind.data.frame(Year = 2020, weight_mean)
colnames(weight_mean)[3] <- "Weight"
weight <- rbind.data.frame(weight, weight_mean)

# Plot
weight$Age <- as.numeric(weight$Age)
weight$Year <- as.numeric(weight$Year)
weight_plot <- ggplot(weight, aes(x = Year, y = Age, fill = Weight)) +
  geom_tile() +
  scale_y_continuous(expand = c(0, 0), 
                     breaks=c(1, 5, 10, 15), 
                     labels = c(1, 5, 10, "15+")) + 
  scale_x_continuous(expand = c(0, 0), 
                     breaks=c(1980, 1990, 2000, 2010, 2020),
                     labels = c(1980, 1990, 2000, 2010, "Average")) +
  geom_vline(xintercept = 2019.5, color = "lightgray") +
  scale_fill_viridis(name = "Weight (kg)", limits = c(0, NA)) +
  coord_equal() +
  ylab("Age") + xlab("Year") +
  theme(panel.border = element_rect(colour = NA, fill = NA))
weight_plot

# Combine sensitivity popdy data and plot together ----------------------------
diet_popdy <- read.csv("sensitivities/diet_popdy.csv")
time_popdy <- read.csv("sensitivities/time_popdy.csv")

sens_popdy <- rbind.data.frame(cbind.data.frame(diet_popdy, 
                                                Sensitivity = "Diet Proportion"),
                               cbind.data.frame(time_popdy, 
                                                Sensitivity = "Time Period")) %>%
  filter(year <= 2022) %>%
  mutate(model = factor(model, 
                        levels = c("Max 0.05", "Max 0.1", "Max 0.5", "Max 0.75",
                                   "Base Cannibalism",
                                   "High (1988-1999)", "Low (2005-2019)"))) %>%
  mutate(variable = factor(variable, 
                           levels = c("Spawning Biomass (Mt)", 
                                      "Total Biomass (Mt)", 
                                      "Recruitment (millions)"))) %>%
  ggplot(., aes(x=year, y=value, color = model, fill = model)) +
  geom_line(aes(linetype = model)) +
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed", "solid", "solid")) +
  geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
  scale_color_manual(values = c("#78D6AEFF", "#38AAACFF", "#357BA2FF", "#40498EFF", 
                                "gray20", 
                                "#781C6DFF", "#ED6925FF")) +
  scale_fill_manual(values = c("#78D6AEFF", "#38AAACFF", "#357BA2FF", "#40498EFF", 
                               "gray20", 
                               "#781C6DFF", "#ED6925FF")) +
  # scale_color_viridis(discrete = TRUE, direction = -1) +  
  # scale_fill_viridis(discrete = TRUE, direction = -1) +  
  geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  ylim(0, NA) +
  xlim(1980, 2022) +
  ylab(" ") + xlab("Year") +
  labs(color = "Model", fill = "Model", linetype = "Model") +
  facet_grid(variable ~ Sensitivity, scales = "free_y", switch = "y") +
  theme(strip.background = element_blank(), strip.placement = "outside") 
sens_popdy


# Predicted annual diet composition? ------------------------------------------
pred_diet <- as.data.frame(as.table(ms_priorM1$model$quantities$diet_prop_weight_hat)) %>%
  filter(Var1 == "A" & Var2 == "A")  # filter to predator/prey sex = A

# Update labeling for ages and years, rename columns
levels(pred_diet$Var3) <- c(1:20)  
levels(pred_diet$Var4) <- c(1:20)
levels(pred_diet$Var5) <- c(all_yrs)

pred_diet <- pred_diet[, -(1:2)]
colnames(pred_diet) <- c("Pred_age", "Prey_age", "Year", "proportion")
pred_diet <- pred_diet %>% 
  mutate(Year = as.numeric(as.character(Year))) %>%
  filter(Year < 2020) %>%  # remove projection
  mutate(Type = "Predicted")  # add label column

# Read in annual diet proportions 
# TODO: double check that these are ok!
annual_diet <- read.csv("data/diet/diet_for_CEATTLE_yearly.csv")[, c("Pred_age", "Prey_age", "Year",
                                                                     "Stomach_proportion_by_weight")]
annual_diet$Type <- "Observed"
colnames(annual_diet)[4] <- "proportion"
 
diet_compared <- rbind.data.frame(pred_diet, annual_diet) %>%
  mutate(Prey_age = as.numeric(Prey_age),
         Pred_age = as.numeric(Pred_age)) %>%
  filter(Prey_age < 6) %>%  # remove prey above age 5 (no predation on those ages)
  filter(Pred_age < 16) %>%  # remove predators above age 15 (as it's a plus group)
  mutate(Prey_age = factor(Prey_age),
         Pred_age = factor(Pred_age)) %>%
  ggplot(., aes(x=Year, y=proportion, fill = Prey_age)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    ylab("Diet Proportion by Weight") + labs(fill = "Prey Age") +
    facet_grid(Pred_age ~ Type)
diet_compared

# Read in number of samples & plot
sample_size <- read.csv("data/diet/diet_for_CEATTLE_yearly.csv")[, c("Pred_age", "Prey_age", "Year",
                                                                     "Sample_size")] %>%
  mutate(Type = rep("Sample Size")) %>%
  mutate(Prey_age = as.numeric(Prey_age),
         Pred_age = as.numeric(Pred_age)) %>%
  filter(Prey_age < 6) %>%  # remove prey above age 5 (no predation on those ages)
  filter(Pred_age < 16) %>%  # remove predators above age 15 (as it's a plus group)
  mutate(Prey_age = factor(Prey_age),
         Pred_age = factor(Pred_age)) %>%
  ggplot(., aes(x=Year, y=Sample_size, fill = Prey_age)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  ylab("Number of Prey Items") + labs(fill = "Prey Age") +
  facet_grid(Pred_age ~ Type) +
  theme(legend.position = "none")
sample_size

predict_diet_plot <- cowplot::plot_grid(sample_size, diet_compared,
                                        rel_widths = c(1, 2.1))
predict_diet_plot


### Save plots (when not experimenting) ---------------------------------------
# ggsave(filename="plots/CEATTLE/cannibalism/popdyn_M1prior.png", plots$popdy, width=170, height=180, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/biomass_ratio.png", plots$ratio, width=150, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/nbyage.png", plots$nbyage, width=160, height=120, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/nbyage_anomaly.png", plots$nbyage_anomaly, width=170, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/biom_n_byage.png", biom_n_byage, width=170, height=170, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/survey_biomass.png", plots$survey, width=200, height=120, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/suitability.png", plots$suit, width=150, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/biomass_byage.png", plots$biombyage, width=160, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/biomass_consumed.png", plots$b_consumed, width=140, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/realized_consumption.png", plots$yearly_b, width=140, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/M.png", ms_prior_mort[[1]], width = 160, height = 70, units = "mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/F.png", eff_plot, width=150, height=80, units="mm", dpi=300)
# ggsave(filename="plots/weight-at-age.png", weight_plot, width = 180, height = 70, units = "mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/Testing/ALL_sens_popdy.png", sens_popdy, width = 170, height = 120, units = "mm", dpi = 300)
# ggsave(filename="plots/CEATTLE/predicted_diet_comparison.png", predict_diet_plot, width = 230, height = 280, units = "mm", dpi = 300)
