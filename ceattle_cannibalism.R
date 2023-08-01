#' Plot results of CEATTLE runs (from run_ceattle.R), comparing the single-
#' species model and the multispecies model (cannibalism), which uses 
#' proportion of cannibalism calculated from diet database going back to 1988. 
#' Different specifications of M1 (fixed, estimated, or with a prior) are also 
#' explored.

# devtools::install_github("grantdadams/Rceattle", ref = "dev")
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggview)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

### Load in models - Rdata objects --------------------------------------------
load("models/ss_fixM1.Rdata")
load("models/ss_estM1.Rdata")
load("models/ss_priorM1.Rdata")
load("models/ms_estM1.Rdata")
load("models/ms_fixM1.Rdata")
load("models/ms_priorM1.Rdata")

# Compare fits
model_fits <- rbind(cbind(model = "SS est M1", ss_estM1$fit),
                    cbind(model = "SS fix M1", ss_fixM1$fit),
                    cbind(model = "SS prior M1", ss_priorM1$fit),
                    cbind(model = "MS est M1", ms_estM1$fit),
                    cbind(model = "MS fix M1", ms_fixM1$fit),
                    cbind(model = "MS prior M1", ms_priorM1$fit))
# model_summary <- cbind(ss_estM1$summary,
#                        ss_fixM1$summary[, 3],
#                        ss_priorM1$summary[, 3],
#                        ms_estM1$summary[, 3],
#                        ms_fixM1$summary[, 3],
#                        ms_priorM1$summary[, 3])[, -(1:2)]
# colnames(model_summary) <- c("SS est M1", "SS fix M1", "SS prior M1",
#                              "MS est M1", "MS fix M1", "MS prior M1")

# New JNLL component tables (until models get re-run with this code included)
comp_out <- function(run) {
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
  return(comp)
}
model_summary <- cbind(comp_out(ss_estM1$model),
                       comp_out(ms_estM1$model)[, 2])
colnames(model_summary) <- c("component", "SS model", "MS model")

### Plot multi-species vs. single-species vs. assessment ----------------------
start_yr <- ms_estM1$model$data_list$styr
end_yr <- 2022
years <- start_yr:end_yr
hind_end <- 2019
assess_yr = "2020"

# Helper function for extracting -by-age data from CEATTLE
extract_byage <- function(result, name, type) {
  df <- as.data.frame(as.table(result))
  
  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(years)
  colnames(df) <- c("age", "year", type)
  
  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"
  
  return(df)
}

plot_models <- function(ms_run, ss_run, save_data = FALSE) {
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
  
  # Get recruitment from SS3 files, offset by 1 year
  ss3_R <- read.table(paste0("data/assessment/", assess_yr, "/recruitment.txt"))[-(1:((start_yr-1)-start)), ]
  ss3_R <- ss3_R[-nrow(ss3_R), ]
  
  # Find mean difference between the model runs -------------------------------
  rel_change <- function(df1, df2, title) {
    mean_out <- mean((df1 / 1000000) - (df2 / 1000000))
    SEM <- sd((df1 / 1000000) - (df2 / 1000000)) / sqrt(length(range))
    percent <- mean(((df1 - df2) / df2) * 100) 
    return(c(title, mean_out, SEM, percent))
  }
  
  n_row <- length(start_yr:end_yr)
  rechange_all <- rbind(rel_change(biomass[1:n_row, 3],
                                   nodiet_biomass[1:n_row, 3],
                                   "cannibalism - no diet, SSB"),
                        rel_change(biomass[(n_row+1):(2*(n_row)), 3],
                                   nodiet_biomass[(n_row+1):(2*(n_row)), 3],
                                   "cannibalism - no diet, Total"),
                        rel_change(recruitment, nodiet_R,
                                   "cannibalism - no diet, R"),
                        rel_change(nodiet_biomass[1:n_row, 3],
                                   ss3_biom[1:n_row, 3],
                                   "no diet - SS3, SSB"),
                        rel_change(nodiet_biomass[(n_row+1):(2*(n_row)), 3],
                                   ss3_biom[(n_row+1):(2*(n_row)), 3],
                                   "no diet - SS3, Total"),
                        rel_change(nodiet_R,
                                   ss3_R[1:(n_row), 2],
                                   "no diet - SS3, R"))

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
                               variable = rep("Assessment"), 
                               value = ss3_R[, 2],
                               error = ss3_R[, 3]))
  R_all <- rbind(cbind(R, error = c(ms_run$sdrep$sd[which(names(ms_run$sdrep$value) == "R")], 
                                    ss_run$sdrep$sd[which(names(ss_run$sdrep$value) == "R")])), 
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
  
  # Plot ratio of SSB:Biomass to look for skewness in age composition
  ratio <- as.data.frame(cbind(year = years,
                               assessment = ss3_ssb[, 1]/ss3_biomass[, 1],
                               cannibalism = biomass[1:length(years), 3] / 
                                 biomass[(length(years)+1):(length(years)*2), 3],
                               single_species = nodiet_biom[1:length(years), 3] / 
                                 nodiet_biom[(length(years)+1):(length(years)*2), 3]))
  
  colnames(ratio) <- c("year", "Assessment", "CEATTLE - cannibalism", 
                       "CEATTLE - single-species")
  ratio2 <- melt(ratio, id.vars = "year", variable.name = "model")
  ratio2$model <- factor(ratio2$model, levels = c("Assessment", 
                                                  "CEATTLE - single-species", 
                                                  "CEATTLE - cannibalism"))
  
  ratio_plot <- ggplot(ratio2, aes(x=year, y=value, color=model)) +
    geom_line() +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    geom_vline(xintercept = hind_end, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylab("SSB/Biomass")
  
  # Plot the difference between model runs
  ceattle_intrasp <- all_popdy %>% filter(model == "CEATTLE - cannibalism")
  ceattle_nodiet <- all_popdy %>% filter(model == "CEATTLE - single-species")
  assessment <- all_popdy %>% filter(model == "Assessment")

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
    theme(strip.background = element_blank(), strip.placement = "outside")
  
  # Plot numbers-at-age -------------------------------------------------------
  nbyage <- extract_byage(ms_run$quantities$NByage, "CEATTLE - cannibalism", "numbers")
  
  # # Read in data from no diet CEATTLE run
  # nbyage_nodiet <- extract_nbyage(ss_run, "CEATTLE - single species")
  
  # Read in data from SS3 & average beginning & middle of the year
  nbyage_ss3_all <- read.csv(paste0("data/assessment/", assess_yr, "/nbyage.csv")) %>%
    filter(Yr >= start_yr & Yr <= end_yr)
  colnames(nbyage_ss3_all) <- c("year", "timing", c(0:20))
  
  nbyage_ss3_wide <- nbyage_ss3_all %>%
    group_by(year) %>%
    summarize_at(vars("0":"20"), mean)
  
  nbyage_ss3 <- melt(nbyage_ss3_wide[, -2], id.vars = "year")
  nbyage_ss3 <- cbind(nbyage_ss3, 
                      rep("Stock Assessment", length(nbyage_ss3$year)))
  colnames(nbyage_ss3)[2:4] <- c("age", "numbers", "model")
  
  # Combine with nbyage from intrasp run
  nbyage_all <- rbind(nbyage, nbyage_ss3)
  
  # Find mean numbers-at-age for each year
  nbyage_all$age <- as.numeric(nbyage_all$age)
  nbyage_mean <- nbyage_all %>%
    group_by(year, model) %>%
    summarize(mean = weighted.mean(age, numbers))
  nbyage_mean$year <- as.numeric(nbyage_mean$year)
  
  mean_nbyage_plot <- ggplot(nbyage_mean, aes(x = year, y = mean)) +
    geom_line(aes(color = model)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) 
  
  # Set 15 as accumulation age
  nbyage_all$age[as.numeric(nbyage_all$age) > 15] <- 15
  
  # Reduce down to millions for plotting
  nbyage_all$numbers <- nbyage_all$numbers / 1000000
  
  # Plot yearly nbyage
  nbyage_plot <- ggplot(nbyage_all, aes(x=year, y=age)) +
    geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
    scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
    scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
    scale_y_continuous(breaks = seq(1, 15, 2), labels = c(seq(1, 13, 2), "15+")) +
    scale_x_discrete(breaks = seq(start_yr, end_yr, 3)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    xlab(" ") + ylab("Age") + 
    labs(fill="millions (n)", size="millions (n)", color="millions (n)") +
    facet_wrap(~model, ncol=1)
  
  # Plot comparison to survey index -------------------------------------------
  init_surv <- ms_estM1$model$data_list$srv_biom %>% filter(Year > 1)  # input survey biomass
  
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
  survey <- cbind(survey, model = "Assessment")
  
  survey_all <- rbind(intrasp_srv, nodiet_srv, survey)
  survey_all$model <- factor(survey_all$model, levels = c("Assessment", "CEATTLE - single-species", "CEATTLE - cannibalism"))
  
  survey_plot <- ggplot() +
    geom_vline(xintercept = hind_end, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    geom_pointrange(data = init_surv, 
                    aes(x = Year, y = Observation,
                        ymin = exp(log(Observation) - 1.96*Log_sd),
                        ymax = exp(log(Observation) + 1.96*Log_sd)),
                    fatten = 5) +
    geom_line(data = survey_all, aes(x = year, y = biomass, color = model), alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    ylim(0, NA) +
    xlab("year") + ylab("Index of Abundance")
  
# Suitability ---------------------------------------------------------------
  suitability <- as.data.frame(as.table(ms_run$quantities$suit_main)) %>%
    filter(Var1 == "A" & Var2 == "A")
  
  suitability <- suitability[, 3:6]
  suitability$Var3 <- factor(suitability$Var3, labels = c(1:15))
  suitability$Var4 <- as.integer(factor(suitability$Var4, labels = c(1:15)))
  suitability$Var5 <- factor(suitability$Var5, labels = years)
  
  colnames(suitability) <- c("pred_age", "prey_age", "year", "value")
  suitability <- suitability %>% filter(prey_age < 6)
  suitability$prey_age <- factor(suitability$prey_age)
  suitability <- suitability %>% 
    group_by(pred_age, prey_age) %>%
    summarize(value = mean(value))
  
  suit_plot <- ggplot(suitability, aes(x = pred_age, y = value, fill = prey_age)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, name = "prey age") +
    xlab("predator age") + ylab("suitability") 
  
  # Biomass-at-age for each model run -----------------------------------------
  biombyage <- extract_byage(ms_run$quantities$biomassByage, 
                             "CEATTLE - cannibalism", "biomass")
  
  # Set 15 as accumulation age
  biombyage$age[as.numeric(biombyage$age) > 15] <- 15
  biombyage$age <- as.integer(biombyage$age)
  
  # Plot yearly biomass by age
  biombyage_plot <- ggplot(biombyage, aes(x=year, y=age)) +
    geom_point(aes(size = biomass, color = biomass, fill = biomass)) +
    scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
    scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
    scale_y_continuous(breaks = seq(1, 15, 2), labels = c(seq(1, 13, 2), "15+")) +
    scale_x_discrete(breaks = seq(start_yr, end_yr, 3)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    xlab(" ") + ylab("Age") 
  
  ### Plot realized consumption -------------------------------------------------
  # Extract biomass consumed as prey
  b_consumed <- extract_byage(ms_run$quantities$B_eaten_as_prey, 
                              "CEATTLE - cannibalism", "biomass")[, -4]
  b_consumed$age <- as.integer(b_consumed$age)
  
  # Filter to only ages below 6 (the ages of hake consumed)
  b_consumed <- b_consumed %>% filter(age < 6)
  b_consumed$age <- as.factor(b_consumed$age)
  
  # Plot yearly biomass consumed by age
  b_consumed_plot <- ggplot(b_consumed, aes(x=year, y=(biomass / 1000000), fill = age)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    scale_x_discrete(breaks = seq(start_yr, end_yr, 3)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylab("Biomass consumed (Mt)")
  
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
    geom_vline(xintercept = hind_end, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylab("biomass of prey / SSB")
  
  if(save_data == TRUE) {
    write.csv(biomass, paste0("data/CEATTLE/", assess_yr, "/ceattle_intrasp_biomass.csv"), row.names = FALSE)
    write.csv(nodiet_biomass, paste0("data/CEATTLE/", assess_yr, "/ceattle_nodiet_biomass.csv"), row.names = FALSE)
    write.csv(recruitment, paste0("data/CEATTLE/", assess_yr, "/ceattle_intrasp_R.csv"), row.names = FALSE)
    write.csv(nbyage, paste0("data/CEATTLE/", assess_yr, "/ceattle_intrasp_nbyage.csv"), row.names = FALSE)
  }
  
  return(list(relative_change = rechange_all, popdy = popdy_plot, ratio = ratio_plot, 
              pop_diff = diff_plot, nbyage = nbyage_plot, survey = survey_plot, 
              suit = suit_plot, biombyage = biombyage_plot, 
              b_consumed = b_consumed_plot, yearly_b = yearly_b_plot))
}

plots <- plot_models(ms_estM1$model, ss_estM1$model)
relative_change <- plots$relative_change
plots$popdy
plots$ratio
plots$pop_diff
plots$nbyage
plots$survey
plots$suit
plots$biombyage
plots$b_consumed
plots$yearly_b

# Plot with fixed M1
plots_M1fixed <- plot_models(ms_fixM1$model, ss_fixM1$model)
relative_change_M1fix <- plots_M1fixed$relative_change
plots_M1fixed$popdy


# Plot with a prior on M1
plots_M1prior <- plot_models(ms_priorM1$model, ss_priorM1$model)
relative_change_M1prior <- plots_M1prior$relative_change
plots_M1prior$popdy

### Compare and plot natural mortality (M1 + M2) ------------------------------
mortality <- function(run, type) {
  M1 <- run$quantities$M1[1, 1, 1:15]
  if(type == "single-species") {
    return(M1)
  }
  
  if(type == "multi-species") {
    M2 <- extract_byage(run$quantities$M2, "multispecies", "M2")
    total_mortality <- M2 %>%
      mutate(M1_M2 = M2 + rep(M1, length(years)))
    total_mortality$age <- as.integer(total_mortality$age)
    total_mortality$year <- as.integer(as.character(total_mortality$year))
    
    mortality_plot <- ggplot(total_mortality, aes(y = age, x = year, zmin = 0, zmax = 1.6)) + 
      geom_tile(aes(fill = M1_M2)) +
      scale_y_continuous(expand = c(0, 0), breaks=c(1, 3, 5, 7, 9, 11, 13, 15)) + 
      scale_x_continuous(expand = c(0, 0), breaks=c(1990, 1995, 2000, 2005, 2010, 2015, 2020)) + 
      scale_fill_viridis(name = "M1 + M2", limits = c(0, 2.1), breaks = c(0.21, 1, 2)) +
      geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
      coord_equal() +
      ylab("Age") + xlab("Year") +
      theme(panel.border = element_rect(colour = NA, fill = NA))
    
    # Min, max, mean natural mortality by age, for ages 1-5
    M_byage <- total_mortality %>%
      filter(age < 6) %>%
      group_by(age) %>%
      summarize(min = min(M1_M2), max = max(M1_M2), mean = mean(M1_M2))
    
    return(list(mortality_plot, M_byage, M1))
  }
}

# Cannibalism with estimated M1
ms_mort <- mortality(ms_estM1$model, type = "multi-species")
ms_mort[[1]]
ms_mort[[2]]
ms_M1 <- ms_mort[[3]]

# Cannibalism with fixed M1
ms_fixed_mort <- mortality(ms_fixM1$model, type = "multi-species")
ms_fixed_mort[[1]]
ms_fixed_mort[[2]]

# Cannibalism with prior on M1
ms_prior_mort <- mortality(ms_priorM1$model, type = "multi-species")
ms_prior_mort[[1]]
ms_prior_mort[[2]]
ms_prior_M1 <- ms_prior_mort[[3]]

# Single-species with estimated M1
ss_M1 <- mortality(ss_estM1$model, type = "single-species")
ss_prior_M1 <- mortality(ss_priorM1$model, type = "single-species")

### Reference points ----------------------------------------------------------
# Catch / Fishing effort
eff <- rbind(data.frame(year = years,
                        eff = t(ms_fixM1$model$quantities$F_spp),
                        model = "cannibalism"),
             data.frame(year = years,
                        eff = t(ss_fixM1$model$quantities$F_spp),
                        model = "single-species"))

eff_plot <- ggplot(eff, aes(x=year, y=eff, color=model)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  xlab("year") + ylab("F")
eff_plot


ref_points <- cbind.data.frame(metric = c("R0", "SPR at B40"),
                               SS = c((ss_estM1$model$quantities$R0 / 100000),
                                      ss_estM1$model$quantities$SPRlimit),
                               MS = c((ms_estM1$model$quantities$R0 / 100000),
                                      ms_estM1$model$quantities$SPRlimit))

quantile((ss_estM1$model$quantities$DynamicSB0 / 1000000), probs = c(0.025, 0.5, 0.975))
quantile((ms_estM1$model$quantities$DynamicSB0 / 1000000), probs = c(0.025, 0.5, 0.975))

# ms_estM1$model$quantities$Flimit # F that gives SPR40%
# ss_estM1$model$quantities$Flimit # F that gives SPR40%
# ms_estM1$model$quantities$Ftarget
# ss_estM1$model$quantities$Ftarget


# ms_run_Fspr <- Rceattle::fit_mod(data_list = ms_estM1$model$data_list,
#                                  inits =  ms_estM1$model$estimated_params, # Initial parameters from ss_run
#                                  estimateMode = 2, # Run projection only
#                                  HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
#                                                            FsprLimit = 0.4, # F40%
#                                                            Ptarget = 0.4, # Target is 40% B0
#                                                            Plimit = 0.1, # No fishing when SB<SB10
#                                                            Pstar = 0.45,
#                                                            Sigma = 0.5),
#                                  msmMode = 1, # Single species mode
#                                  verbose = 1)
# 
# ms_run_Fspr$quantities$DynamicB0 #B0
# 
# ms_run_Fspr$quantities$DynamicSB0 #SB0
# 
# ms_run_Fspr$quantities$Flimit #F that gives you SPR40%
# 
# ms_run_Fspr$quantities$SPRlimit #SPR40%

# Relative SSB
relativeSSB_estM1 <- data.frame(t(ms_estM1$model$quantities$biomassSSB / ms_estM1$model$quantities$DynamicSB0))
relativeSSB_estM1$year <- rownames(relativeSSB_estM1)
rownames(relativeSSB_estM1) <- NULL
relativeSSB_estM1$model <- "estimated M1"
quantile(relativeSSB_estM1$Hake, probs = c(0.025, 0.5, 0.975))

relativeSSB_fixM1 <- data.frame(t(ms_fixM1$model$quantities$biomassSSB / ms_fixM1$model$quantities$DynamicSB0))
relativeSSB_fixM1$year <- rownames(relativeSSB_fixM1)
rownames(relativeSSB_fixM1) <- NULL
relativeSSB_fixM1$model <- "fixed M1"

relativeSSB_priorM1 <- data.frame(t(ms_priorM1$model$quantities$biomassSSB / ms_priorM1$model$quantities$DynamicSB0))
relativeSSB_priorM1$year <- rownames(relativeSSB_priorM1)
rownames(relativeSSB_priorM1) <- NULL
relativeSSB_priorM1$model <- "prior M1"

relativeSSB <- rbind(relativeSSB_estM1, relativeSSB_fixM1, relativeSSB_priorM1)
colnames(relativeSSB)[1] <- "relativeSSB"
relativeSSB$year <- as.numeric(relativeSSB$year)
relativeSSB$model <- factor(relativeSSB$model)
ggplot(relativeSSB, aes(x = year, y = relativeSSB, color = model)) +
  geom_line() +
  geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  ylab("Relative SSB")


### Save plots (when not experimenting) ---------------------------------------
ggsave(filename="plots/CEATTLE/cannibalism/popdyn.png", plots$popdy, width=140, height=150, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/biomass_ratio.png", plots$ratio, width=150, height=80, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/nbyage.png", plots$nbyage, width=160, height=120, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/survey_biomass.png", plots$survey, width=200, height=120, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/suitability.png", plots$suit, width=150, height=80, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/biomass_byage.png", plots$biombyage, width=160, height=80, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/biomass_consumed.png", plots$b_consumed, width=140, height=80, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/realized_consumption.png", plots$yearly_b, width=140, height=80, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/M.png", ms_mort[[1]], width = 160, height = 70, units = "mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/popdyn_M1fixed.png", plots_M1fixed$popdy, width=140, height=150, units="mm", dpi=300)
ggsave(filename="plots/CEATTLE/cannibalism/popdyn_M1prior.png", plots_M1prior$popdy, width=140, height=150, units="mm", dpi=300)
