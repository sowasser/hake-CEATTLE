# Run CEATTLE with intraspecies-predation proportions calculated from diet 
# database going back to 1980.

# devtools::install_github("grantdadams/Rceattle", ref = "dev")
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggview)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Read in CEATTLE data from the excel file
hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230324.xlsx")

# Run and fit the CEATTLE model -----------------------------------------------
run_CEATTLE <- function(data, M1, init, msm) {
  data$est_M1 <- M1  # Set M1 to fixed (0), estimated (1), age-varying estimate (3)
  # data$endyr <- 2019
  prior = FALSE
  run <- Rceattle::fit_mod(data_list = data,
                           inits = init,
                           file = NULL, # Don't save
                           # debug = 1, # 1 = estimate, 0 = don't estimate
                           msmMode = msm, # Single-species mode - no predation mortality
                           # M1Fun = Rceattle::build_M1(M1_model = M1, 
                           #                            updateM1 = TRUE,
                           #                            M1_use_prior = prior,
                           #                            M1_prior_mean = 0.2,
                           #                            M1_prior_sd = .1),
                           # proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
                           estimateMode = 0,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
                           HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                                                     FsprLimit = 0.4, # F40%
                                                     Ptarget = 0.4, # Target is 40% B0
                                                     Plimit = 0.1, # No fishing when SB<SB10
                                                     Pstar = 0.45,
                                                     Sigma = 0.5),
                           phase = "default",
                           # initMode = 1,
                           verbose = 1)
  
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  jnll_summary <- as.data.frame(run$quantities$jnll_comp)
  jnll_summary$sum <- rowSums(run$quantities$jnll_comp)
  return(list(run, fit, jnll_summary))
}

# No diet (single-species run)
nodiet <- run_CEATTLE(data = hake_intrasp, M1 = 1, init = NULL, msm = 0)
nodiet[[2]]  # check convergence

# Run with cannibalism, estimated M1
# This uses estimated parameters from the single-species run to help convergence
intrasp <-  run_CEATTLE(data = hake_intrasp, M1 = 1, init = nodiet[[1]]$estimated_params, msm = 1)
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

# Run with cannibalism, fixed M1 & time-varying M1
intrasp_M1fixed <- run_CEATTLE(data = hake_intrasp, M1 = 0, init = NULL, msm = 1)
intrasp_M1aged <- run_CEATTLE(data = hake_intrasp, M1 = 3, init = NULL, msm = 1)

# Compare fits
intrasp_fit <- rbind(cbind(model = "est M1", intrasp[[2]]),
                     cbind(model = "fix M1", intrasp_M1fixed[[2]]),
                     cbind(model = "age M1", intrasp_M1aged[[2]]))
intrasp_summary <- cbind(intrasp[[3]], 
                         intrasp_M1fixed[[3]][, 3], 
                         intrasp_M1aged[[3]][, 3])[, -(1:2)]
colnames(intrasp_summary) <- c("est M1", "fix M1", "age M1")

# # Trying age-varying M1 with 3 age blocks: 1, 2, and 3+
# map <- intrasp_run$map
# orig <- intrasp_run$map$mapList$ln_M1
# orig[orig > 2] <- 3
# map$mapList$ln_M1 <- orig
# map$mapFactor$ln_M1 = as.factor(map$mapList$ln_M1)
# 
# intrasp_run2 <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                   inits = NULL,
#                                   file = NULL, # Don't save
#                                   map = map,
#                                   msmMode = 1, # Single-species mode - no predation mortality
#                                   phase = "default")

# Run single-species model with fixed & time-varying M1
nodiet_M1fixed <- run_CEATTLE(data = hake_intrasp, M1 = 0, init = NULL, msm = 0)
nodiet_M1aged <- run_CEATTLE(data = hake_intrasp, M1 = 3, init = NULL, msm = 0)

nodiet_fit <- rbind(cbind(model = "est M1", nodiet[[2]]),
                    cbind(model = "fix M1", nodiet_M1fixed[[2]]),
                    cbind(model = "age M1", nodiet_M1aged[[2]]))

nodiet_summary <- cbind(nodiet[[3]], 
                        nodiet_M1fixed[[3]][, 3],
                        nodiet_M1aged[[3]][, 3])[, -(1:2)]
colnames(nodiet_summary) <- c("est M1", "fix M1", "age M1")



### Plot multi-species vs. single-species vs. assessment ----------------------
start_yr <- intrasp[[1]]$data_list$styr
end_yr <- 2022
years <- start_yr:end_yr

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

plot_models <- function(ms_run, ss_run, assess_yr = "2020", hind_end = 2019, save_data = FALSE) {
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
  
  # Find mean difference between the model runs -------------------------------
  mean_SEM <- function(df1, df2, title) {
    mean_out <- mean((df1 / 1000000) - (df2 / 1000000))
    SEM <- sd((df1 / 1000000) - (df2 / 1000000)) / sqrt(length(range))
    return(c(title, mean_out, SEM))
  }
  
  n_row <- length(start_yr:end_yr)
  mean_SEM_all <- rbind(mean_SEM(biomass[1:n_row, 3],
                                 nodiet_biomass[1:n_row, 3],
                                 "cannibalism - no diet, SSB"),
                        mean_SEM(biomass[(n_row+1):(2*(n_row)), 3],
                                 nodiet_biomass[(n_row+1):(2*(n_row)), 3],
                                 "cannibalism - no diet, Total"),
                        mean_SEM(recruitment, nodiet_R,
                                 "cannibalism - no diet, R"),
                        mean_SEM(nodiet_biomass[1:n_row, 3],
                                 ss3_biom[1:n_row, 3],
                                 "no diet - SS3, SSB"),
                        mean_SEM(nodiet_biomass[(n_row+1):(2*(n_row)), 3],
                                 ss3_biom[(n_row+1):(2*(n_row)), 3],
                                 "no diet - SS3, Total"),
                        mean_SEM(nodiet_R,
                                 ss3_R[1:(n_row), 2],
                                 "no diet - SS3, R"))

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
  
  # TODO: fix this bit!
  # # Plot the difference between model runs
  # ceattle_intrasp <- all_popdy %>% filter(model == "CEATTLE - cannibalism")
  # ceattle_nodiet <- all_popdy %>% filter(model == "CEATTLE - single-species")
  # assessment <- all_popdy %>% filter(model == "Assessment")
  # 
  # diff_intrasp <- rbind(cbind.data.frame(year = years,
  #                                        type = ceattle_intrasp$type,
  #                                        difference = ceattle_intrasp$value - ceattle_nodiet$value,
  #                                        models = "cannibalism - single-species"),
  #                       cbind.data.frame(year = years,
  #                                        type = ceattle_nodiet$type,
  #                                        difference = ceattle_nodiet$value - assessment$value,
  #                                        models = "single-species - assessment"))
  # 
  # diff_plot <- ggplot(diff_intrasp, aes(x=year, y=difference, color=models)) +
  #   geom_line() +
  #   scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  #   facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
  #   theme(strip.background = element_blank(), strip.placement = "outside")
  
  
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
  survey_biom <- function(run, name) {
    srv <- data.frame(year = 1995:hind_end,
                      biomass = run$quantities$srv_bio_hat,
                      log_sd = run$quantities$srv_log_sd_hat,
                      model = rep(name, length(1995:hind_end)))
    return(srv)
  }
  
  intrasp_srv <- survey_biom(ms_run, "CEATTLE - cannibalism")
  nodiet_srv <- survey_biom(ss_run, "CEATTLE - no diet")
  
  survey <- read.csv(paste0("data/assessment/", assess_yr, "/survey_data.csv"))
  survey <- cbind(survey, model = rep("Stock Synthesis", length(survey$year)))
  
  survey_all <- rbind(intrasp_srv, nodiet_srv, survey)
  
  survey_plot <- ggplot(survey_all, aes(x=year, y=biomass, color=model)) +
    geom_line(alpha = 0.3) +
    geom_point() +
    # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    geom_vline(xintercept = hind_end, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    xlab("year") + ylab("survey biomass")
  
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
  
  return(list(mean_SEM_all, popdy_plot, ratio_plot, nbyage_plot, 
              survey_plot, suit_plot, biombyage_plot, b_consumed_plot, 
              yearly_b_plot))  # diff_plot
}

plots <- plot_models(intrasp[[1]], nodiet[[1]])
mean_SEM <- plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]
plots[[6]]
plots[[7]]
plots[[8]]
plots[[9]]
# plots[[10]]

# Plot with fixed M1
plots_M1fixed <- plot_models(intrasp_M1fixed[[1]], nodiet_M1fixed[[1]])
plots_M1fixed[[2]]


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
    
    mortality_plot <- ggplot(total_mortality, aes(y = age, x = year, zmin = 0, zmax = 1.5)) + 
      geom_tile(aes(fill = M1_M2)) +
      scale_y_continuous(expand = c(0, 0), breaks=c(1, 3, 5, 7, 9, 11, 13, 15)) + 
      scale_x_continuous(expand = c(0, 0), breaks=c(1990, 1995, 2000, 2005, 2010, 2015, 2020)) + 
      scale_fill_viridis(name = "M1 + M2", limits = c(0, 2.0), breaks = c(0.21, 2.0)) +
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
intrasp_mort <- mortality(intrasp[[1]], type = "multi-species")
intrasp_mort[[1]]
intrasp_mort[[2]]
intrasp_M1 <- intrasp_mort[[3]]

# Cannibalism with fixed M1
intrasp_fixed_mort <- mortality(intrasp_M1fixed[[1]], type = "multi-species")
intrasp_fixed_mort[[1]]
intrasp_fixed_mort[[2]]

# Cannibalism with age-varying M1
intrasp_aged_mort <- mortality(intrasp_M1aged[[1]], type = "multi-species")
intrasp_aged_mort[[1]]
intrasp_aged_mort[[2]]
intrasp_aged_M1 <- intrasp_aged_mort[[3]]

# Single-species with estimated M1
nodiet_M1 <- mortality(nodiet[[1]], type = "single-species")

# Single-species with age-varying M1
nodiet_aged_M1 <- mortality(nodiet_M1aged[[1]], type = "single-species")

### Reference points ----------------------------------------------------------
# TODO: fix this bit! Get specific HCR to work.
# proj_intrasp <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                   inits = intrasp[[1]]$estimated_params, 
#                                   estimateMode = 0,  
#                                   HCR = Rceattle::build_hcr(HCR = 3, 
#                                                             FsprTarget = 0.4, # 0.75 * F40%
#                                                             # FsprLimit = 0.4,
#                                                             Plimit = 0.2))

test <- Rceattle::fit_mod(data_list = hake_intrasp,
                          estimateMode = 0)
Rceattle::plot_recruitment(test, add_ci = TRUE, incl_proj = TRUE)

# Relative SSB
DynamicB0 <- cbind.data.frame(year = years, 
                              B0 = intrasp[[1]]$quantities$DynamicB0[1:length(years)])
SSB <- biomass %>% filter(type == "SSB")
relativeSSB <- cbind.data.frame(year = years, 
                                relativeSSB = SSB$value / 
                                  intrasp[[1]]$quantities$DynamicB0[1:length(years)])
relativeSSB[31,2]  # 2019 value

# Run CEATTLE with an HCR to see reference points
intrasp_Fspr <- Rceattle::fit_mod(data_list = intrasp[[1]]$data_list,
                                  inits = intrasp[[1]]$estimated_params, 
                                  estimateMode = 1,  # hindcast model only
                                  HCR = Rceattle::build_hcr(HCR = 3, 
                                                            FsprTarget = 0.4, # 0.75 * F40%
                                                            # FsprLimit = 0.4,
                                                            Plimit = 0.2))
intrasp_Fspr$quantities$B0
intrasp_Fspr$quantities$SB0
intrasp_Fspr$quantities$Flimit # F that gives SPR40%
intrasp_Fspr$quantities$SPRlimit  # SPR40%, which is the same as B40% b/c no stock-recruit curve
intrasp_Fspr$quantities$Ftarget


### Save plots (when not experimenting) ---------------------------------------
# ggsave(filename="plots/CEATTLE/cannibalism/popdyn.png", plots[[2]], width=140, height=150, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/biomass_ratio.png", plots[[3]], width=150, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/nbyage.png", plots[[4]], width=160, height=120, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/survey_biomass.png", plots[[5]], width=200, height=120, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/suitability.png", plots[[6]], width=150, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/biomass_byage.png", plots[[7]], width=160, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/biomass_consumed.png", plots[[8]], width=140, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/realized_consumption.png", plots[[9]], width=140, height=80, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/M.png", intrasp_mort[[1]], width = 160, height = 70, units = "mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/popdyn_M1fixed.png", plots_M1fixed[[2]], width=140, height=150, units="mm", dpi=300)