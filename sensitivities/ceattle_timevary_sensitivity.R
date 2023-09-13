# CEATTLE for Pacific hake while testing different proportions of cannibalism
# from different time periods, corrected with a Dirichlet multinomial 
# distribution.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Models with M1 prior
load("models/ms_priorM1.Rdata")
# Read in different time period models (specified in run_ceattle.R)
load("models/sensitivity/time-varying/run_90s_prior.Rdata")
load("models/sensitivity/time-varying/run_recent_prior.Rdata")

# Year ranges
all_years <- 1980:2050
no_proj <- 1980:2019

timevary_fits <- rbind(cbind(model = "MS", ms_priorM1$fit),
                       cbind(model = "High (90s)", run_90s_prior$fit),
                       cbind(model = "Low (recent)", run_recent_prior$fit))

timevary_summary <- list(run_90s_prior$summary, run_recent_prior$summary) %>%
  reduce(full_join, by = "component")
colnames(timevary_summary) <- c("component", "High (90s)", "Low (recent)")

### Plot population dynamics --------------------------------------------------
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
  
  # Find mean difference between the model runs
  rel_change <- function(model1, model2, stat, years) {
    df1 <- all_popdy %>% 
      filter(model == model1 & variable == stat & year %in% years)
    df2 <- all_popdy %>%
      filter(model == model2) %>% filter(variable == stat) %>% filter(year %in% years)
    mean_out <- mean((df1$value) - (df2$value))
    SEM <- sd((df1$value) - (df2$value)) / sqrt(length(range))
    percent <- mean(((df1$value - df2$value) / df2$value) * 100) 
    label <- paste(model1, "-", model2, ", ", stat)
    return(c(label, mean_out, SEM, percent))
  }
  
  rechange_all <- rbind(rel_change("High (1988-1999)", "Base Cannibalism Model", "SSB", 1988:1999),
                        rel_change("High (1988-1999)", "Base Cannibalism Model", "Total Biomass", 1988:1999),
                        rel_change("High (1988-1999)", "Base Cannibalism Model", "Recruitment", 1988:1999),
                        rel_change("Low (2005-2019)", "Base Cannibalism Model", "SSB", 2005:2019),
                        rel_change("Low (2005-2019)", "Base Cannibalism Model", "Total Biomass", 2005:2019),
                        rel_change("Low (2005-2019)", "Base Cannibalism Model", "Recruitment", 2005:2019))
  colnames(rechange_all) <- c("model", "mean", "SEM", "percent")
  
  # Plot population dynamics
  all_popdy$variable <- factor(all_popdy$variable, labels = c("SSB (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  all_popdy$model <- factor(all_popdy$model, levels = c("Low (2005-2019)", "High (1988-1999)", "Base Cannibalism Model"))
  
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line() +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    xlim(1980, 2022) +
    ylab(" ") + xlab("Year") +
    labs(color = "Model", fill = "Model") +
    facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")

  return(list(all_popdy, rechange_all, popdy_plot))
}

timing_popdy <- timing_plot_popdy(run_high = run_90s_prior$model, 
                                  run_low = run_recent_prior$model,
                                  ms_model = ms_priorM1$model)
relative_change <- timing_popdy[[2]]
timing_popdy[[3]]

# # Calculate reference points 
# DynamicB0_recent <- c(run_recent_noproj$quantities$DynamicB0[1:length(2005:2019)])
# SSB_recent <-  run_recent_noproj$quantities$biomassSSB * 2
# 
# SSB_recent[15] / DynamicB0_recent[15] # 2019 value


# # Numbers-at-age for each model run -------------------------------------------
# # Read in data from no diet CEATTLE run
extract_byage2 <- function(quantity, name, years) {
  df <- as.data.frame(as.table(quantity))

  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(years)
  colnames(df) <- c("age", "year", "numbers")

  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"
  # years <- as.character(years)
  # df <- df %>% filter(year %in% years)

  return(df)
}
# 
# nbyage_test_all <- rbind(extract_byage2(run_90s_prior$model$quantities$NByage, 
#                                         "High (1988-1999)", 1988:2019),
#                          extract_byage2(ms_priorM1$model$quantities$NByage, 
#                                         "Base Cannibalism Model", all_years),
#                          extract_byage2(run_recent_prior$model$quantities$NByage, 
#                                         "Low (2005-2019)", no_proj) %>% 
#                            filter(year %in% 2005:2019))
# 
# # # Set 15 as accumulation age
# # nbyage_test_all$age[as.numeric(nbyage_test_all$age) > 15] <- 15
# 
# # Plot yearly nbyage
# nbyage_test_all$age <- as.numeric(nbyage_test_all$age)
# nbyage_test_all$model <- factor(nbyage_test_all$model, 
#                                 levels = c("High (1988-1999)", "Base Cannibalism Model", "Low (2005-2019)"))
# nbyage_test_all$year <- factor(nbyage_test_all$year,
#                                levels = c(all_years))
# 
# test_nbyage_plot <- ggplot(nbyage_test_all, aes(x=year, y=age)) +
#   geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
#   scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
#   scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
#   scale_x_discrete(breaks = seq(1980, 2022, 3)) +
#   xlab(" ") + ylab("Age") +
#   theme(legend.position = "none") +
#   facet_wrap(~model, ncol = 1)
# test_nbyage_plot

### Plot mortality ------------------------------------------------------------
extract_M <- function(run, quantity, name) {
  M1 <- run$quantities$M1[1, 1, 1:20]
  M2 <- extract_byage2(quantity, name, all_years)
  total_mortality <- M2 %>%
    mutate(M1_M2 = M2$numbers + rep(M1, length(all_years)))
  total_mortality$age <- as.integer(total_mortality$age)
  total_mortality$year <- as.integer(as.character(total_mortality$year))
  return(total_mortality)
}

M_all <- rbind(extract_M(run_90s_prior$model, run_90s_prior$model$quantities$M2,
                         "High (1988-1999)") %>% filter(year %in% 1980:1999),
               extract_M(run_recent_prior$model, run_recent_prior$model$quantities$M2,
                         "Low (2005-2019)") %>% filter(year %in% 2005:2019),
               extract_M(ms_priorM1$model, ms_priorM1$model$quantities$M2, 
                         "Base Cannibalism Model") %>% filter(year %in% 1980:2022)) 

max(M_all$M1_M2)  # check max M for plotting
timevary_M <- ggplot(M_all, aes(y = age, x = year, zmin = 0, zmax = 1.5)) +
  geom_tile(aes(fill = M1_M2)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_viridis(name = "M1 + M2", limits = c(0, 2.7), breaks = c(0.21, 1, 2, 3)) +
  geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  coord_equal() +
  ylab("Age") + xlab("Year") +
  facet_wrap(~model, ncol = 1)
timevary_M

# Min, max, mean natural mortality by age, for ages 1-5
M_byage <- M_all %>%
  filter(age < 6) %>%
  group_by(age, model) %>%
  summarize(min = min(M1_M2), max = max(M1_M2), mean = mean(M1_M2))

M1_all <- rbind(data.frame(model = "Base Cannibalism Model", 
                           mean = mean(ms_priorM1$model$quantities$M1[1, 1, 1:20])),
                data.frame(model = "High (1988-1999)", 
                           mean = mean(run_90s_prior$model$quantities$M1[1, 1, 1:20])),
                data.frame(model = "Low (2005-2019)", 
                           mean = mean(run_recent_prior$model$quantities$M1[1, 1, 1:20])))



# # Models with estimated M1 ----------------------------------------------------
# load("models/ms_estM1.Rdata")
# # Read in different time period models (specified in run_ceattle.R)
# load("models/sensitivity/time-varying/run_90s.Rdata")
# load("models/sensitivity/time-varying/run_recent.Rdata")
# 
# timing_popdy_estM1 <- timing_plot_popdy(run_high = run_90s$model, 
#                                         run_low = run_recent$model,
#                                         ms_model = ms_estM1$model)
# relative_change_estM1 <- timing_popdy_estM1[[2]]
# timing_popdy_estM1[[3]]


### Save plots (when not experimenting) ---------------------------------------
# ggsave(filename="plots/CEATTLE/cannibalism/Testing/timevarying_popdy.png", timing_popdy[[3]],
#        width=140, height=150, units="mm", dpi=300)
# ggsave(filename = "plots/CEATTLE/cannibalism/Testing/timevarying_nbyage.png", test_nbyage_plot,
#        width=150, height=150, units="mm", dpi=300)
# ggsave(filename = "plots/CEATTLE/cannibalism/Testing/timevarying_M.png",
#        timevary_M, width=140, height = 170, units = "mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/Testing/timevarying_popdy_estM1.png", timing_popdy_estM1[[3]],
#        width=140, height=150, units="mm", dpi=300)
