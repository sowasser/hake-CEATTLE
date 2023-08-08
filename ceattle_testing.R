library(Rceattle)
library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(ggview)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Read in CEATTLE data from the excel file
hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230616.xlsx")

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
                 initMode = 1,
                 # random_rec = TRUE,
                 # random_sel = TRUE,
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

ms_model <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 1, 
                        prior = TRUE, 
                        init = NULL, 
                        msm = 1, 
                        estMode = 0)
ms_model$fit  # check convergence
ms_model$model$quantities$M1

ss_model <- run_CEATTLE(data = hake_intrasp, 
                        M1 = 1, 
                        prior = TRUE, 
                        init = NULL, 
                        msm = 0, 
                        estMode = 0)
ms_model$fit  # check convergence
ms_model$model$quantities$M1


#' Plot results of CEATTLE runs (from run_ceattle.R), comparing the single-
#' species model and the multispecies model (cannibalism), which uses 
#' proportion of cannibalism calculated from diet database going back to 1988. 
#' Different specifications of M1 (fixed, estimated, or with a prior) are also 
#' explored.

# devtools::install_github("grantdadams/Rceattle", ref = "dev")






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
  
  # Difference between both models
  nbyage_diff <- nbyage_ss3
  nbyage_diff$age <- as.numeric(nbyage_diff$age)
  nbyage_diff <- nbyage_diff %>% filter(as.numeric(age) <= 15)
  nbyage_diff$numbers <- (nbyage$numbers - nbyage_diff$numbers) / 1000000
  nbyage_diff$model <- "CEATTLE - assessment"
  nbyage_diff$year <- factor(nbyage_diff$year)
  
  limit <- max(abs(nbyage_diff$numbers)) * c(-1, 1)
  ggplot(nbyage_diff, aes(x=year, y=age)) +
    geom_point(aes(size = numbers, color = numbers)) +
    scale_color_gradientn(colors = pals::ocean.curl(100), limit = limit) +
    scale_y_continuous(breaks = seq(1, 15, 2), labels = c(seq(1, 13, 2), "15+")) +
    scale_x_discrete(breaks = seq(start_yr, end_yr, 3)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    xlab(" ") + ylab("Age") + 
    labs(size="millions (n)", color="millions (n)") 
  
  
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

plots <- plot_models(ms_estM1$model, ss_priorM1$model)
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
