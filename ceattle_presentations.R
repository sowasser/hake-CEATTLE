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
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent_dark.R")
theme_set(theme_sleek_transparent())

### Load in models - Rdata objects --------------------------------------------
load("models/ss_fixM1.Rdata")
load("models/ss_estM1.Rdata")
load("models/ss_priorM1.Rdata")
load("models/ms_estM1.Rdata")
load("models/ms_fixM1.Rdata")
load("models/ms_priorM1.Rdata")

### Plot multi-species vs. single-species vs. assessment ----------------------
start_yr <- ms_estM1$model$data_list$styr
end_yr <- 2022
years <- start_yr:end_yr
hind_end <- 2019
assess_yr = "2020"
max_age <- ms_estM1$model$data_list$nages

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
  
  # Plot popdy ----------------------------------------------------------------
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line(aes(linetype = model)) +
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
  
  # Reduce down to millions for plotting
  nbyage_all$numbers <- nbyage_all$numbers / 1000000
  
  # Plot yearly nbyage
  nbyage_plot <- ggplot(nbyage_all, aes(x=year, y=age)) +
    geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
    scale_fill_viridis(option = "plasma", begin = 0.2) +
    scale_color_viridis(option = "plasma", begin = 0.2) +
    scale_x_discrete(breaks = c(1980, 1990, 2000, 2010, 2020)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    xlab(" ") + ylab("Age") + 
    labs(fill="millions (n)", size="millions (n)", color="millions (n)") +
    facet_wrap(~model, ncol=1)
  
  # # Difference between both models
  nbyage_diff <- cbind.data.frame(year = nbyage$year,
                                  age = as.numeric(nbyage$age),
                                  numbers = (nbyage$numbers - nbyage_ss3$numbers) / 1000000)
  nbyage_diff$model <- "CEATTLE - assessment"
  nbyage_diff$year <- factor(nbyage_diff$year)
  
  limit <- max(abs(nbyage_diff$numbers)) * c(-1, 1)
  nbyage_anomaly <- ggplot(nbyage_diff, aes(x=year, y=age)) +
    geom_point(aes(size = numbers, color = numbers)) +
    scale_color_gradientn(colors = pals::kovesi.diverging_linear_bjr_30_55_c53(100), limit = limit) +
    scale_x_discrete(breaks = c(1980, 1990, 2000, 2010, 2020)) +
    geom_vline(xintercept = as.character(hind_end), linetype = 2, colour = "gray") +  # Add line at end of hindcast
    xlab(" ") + ylab("Age") +
    labs(size="millions (n)", color="millions (n)")
  nbyage_anomaly
  
  return(list(popdy = popdy_plot, nbyage = nbyage_plot, 
              nbyage_anomaly = nbyage_anomaly))
}

plots <- plot_models(ms_priorM1$model, ss_priorM1$model)
plots$popdy
plots$nbyage
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
      mutate(M1_M2 = M2 + rep(M1, length(years)))
    total_mortality$age <- as.integer(total_mortality$age)
    total_mortality$year <- as.integer(as.character(total_mortality$year))
    
    mortality_plot <- ggplot(total_mortality, aes(y = age, x = year, zmin = 0, zmax = 1.6)) + 
      geom_tile(aes(fill = M1_M2)) +
      scale_y_continuous(expand = c(0, 0), breaks=c(1, 5, 10, 15, 20)) + 
      scale_x_continuous(expand = c(0, 0)) + 
      scale_fill_viridis(name = "M1 + M2", limits = c(0, 1.6), 
                         breaks = c(0.21, 1, 1.5),
                         option = "plasma", begin = 0.2) +
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
  scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.2, end = 0.6) +
  ylab("Relative SSB") +
  ylim(0, NA) +
  geom_hline(yintercept = 1, color = "gray") +
  geom_hline(yintercept = 0.4, color = "gray") +
  annotate("text", x = 1985, y = 0.4, label = "Management target", 
           vjust = -0.5, size = 2.5, color = "gray") +
  geom_hline(yintercept = 0.1, color = "gray") +
  annotate("text", x = 1987.5, y = 0.1, label = "Minimum stock size threshold", 
           vjust = -0.5, size = 2.5, color = "gray") 
relativeSSB_plot


### Bioenergetics (temp-dependent consumption) --------------------------------
proxy_params <- read.csv("data/bioenergetics/proxy_bioen_params.csv")

eq_1 <- proxy_params[c(2:4), ]
eq_2 <- proxy_params[c(1, 5:6), ]
eq_2$Tco <- as.numeric(eq_2$Tco)
eq_2$Tcm <- as.numeric(eq_2$Tcm)

### Temperature-dependent consumption -----------------------------------------
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

# Create overall intraspecies predation dataset -------------------------------
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

### Save plots (when not experimenting) ---------------------------------------
ggsave(filename="plots/presentations/popdyn_M1prior.png", plots$popdy, 
       width=200, height=90, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/nbyage.png", plots$nbyage, 
       width=160, height=120, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/nbyage_anomaly.png", plots$nbyage_anomaly, 
       width=200, height=90, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/M.png", ms_prior_mort[[1]], 
       width = 160, height = 70, units = "mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/relative_SSB.png", relativeSSB_plot, 
       width=150, height=80, units="mm", dpi=300, bg = "transparent")
ggsave(filename="plots/presentations/temp_consumption.png", temp_rate,
       bg = "transparent", width=160, height=90, units="mm", dpi=300)
ggsave(filename="plots/presentations/mean_temp_compared.png", mean_temp_compared,
       bg = "transparent", width=180, height=90, units="mm", dpi=300)
ggsave(filename = "plots/presentations/cannibalism_overall.png", diet_plot, 
       bg = "transparent", width=200, height=120, units="mm", dpi=300)
