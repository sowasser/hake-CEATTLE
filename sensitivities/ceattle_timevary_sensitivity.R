# CEATTLE for Pacific hake while testing different proportions of cannibalism
# from different time periods, corrected with a Dirichlet multinomial 
# distribution.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

load("models/ms_estM1.Rdata")
# Read in different time period models (specified in run_ceattle.R)
load("models/sensitivity/time-varying/run_90s.Rdata")
load("models/sensitivity/time-varying/run_recent.Rdata")

### Plot population dynamics --------------------------------------------------
timing_plot_popdy <- function(run_high, run_low) {
  # Pull out SSB & overall biomass from CEATTLE runs
  ceattle_biomass <- function(run, name, years) {
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
  
  test_biom <- rbind(ceattle_biomass(ms_estM1$model, "all years", 1988:2022),
                     ceattle_biomass(run_high, "high (1988-1999)", 1988:1999),
                     ceattle_biomass(run_low, "low (2005-2019)", 2005:2019))
  
  # Put recruitment together
  ceattle_R <- function(run, name, years) {
    R <- c(run$quantities$R)
    error <- c(run$sdrep$sd[which(names(run$sdrep$value) == "R")])
    R_all <- as.data.frame(cbind(year = years, 
                                 variable = rep("Recruitment"),
                                 value = R, 
                                 error = error, 
                                 model = rep(name)))
    return(R_all)
  }
  R_test <- rbind(ceattle_R(ms_estM1$model, "all years", 1988:2022),
                  ceattle_R(run_high, "high (1988-1999)", 1988:1999),
                  ceattle_R(run_low, "low (2005-2019)", 2005:2019))
  
  
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
    mean_out <- mean((df1$value / 1000000) - (df2$value / 1000000))
    SEM <- sd((df1$value / 1000000) - (df2$value / 1000000)) / sqrt(length(range))
    percent <- mean(((df1$value - df2$value) / df2$value) * 100) 
    label <- paste(model1, "-", model2, ", ", stat)
    return(c(label, mean_out, SEM, percent))
  }
  
  rechange_all <- rbind(rel_change("high (1988-1999)", "all years", "SSB", 1988:1999),
                        rel_change("high (1988-1999)", "all years", "Total Biomass", 1988:1999),
                        rel_change("high (1988-1999)", "all years", "Recruitment", 1988:1999),
                        rel_change("low (2005-2019)", "all years", "SSB", 2005:2019),
                        rel_change("low (2005-2019)", "all years", "Total Biomass", 2005:2019),
                        rel_change("low (2005-2019)", "all years", "Recruitment", 2005:2019))
  
  # Plot population dynamics
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
    geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")

  return(list(all_popdy, rechange_all, popdy_plot))
}

timing_popdy <- timing_plot_popdy(run_high = run_90s$model, 
                                  run_low = run_recent$model)
relative_change <- timing_popdy[[2]]
timing_popdy[[3]]

ggsave(filename="plots/CEATTLE/cannibalism/Testing/dirichlet_popdy.png", timing_popdy[[3]], 
       width=140, height=150, units="mm", dpi=300)

# Calculate reference points 
DynamicB0_recent <- c(run_recent$quantities$DynamicB0[1:length(2005:2019)])
SSB_recent <-  run_recent$quantities$biomassSSB * 2

SSB_recent[15] / DynamicB0_recent[15] # 2019 value


# Numbers-at-age for each model run -------------------------------------------
# Read in data from no diet CEATTLE run
extract_nbyage <- function(run, name, years) {
  df <- as.data.frame(as.table(run$quantities$NByage))
  
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

nbyage_test_all <- rbind(extract_nbyage(run_90s$model, "high (1988-1999)", 1988:1999),
                         extract_nbyage(ms_estM1$model, "all years", 1988:2022),
                         extract_nbyage(run_recent$model, "low (2005-2019)", 2005:2019))

# Set 15 as accumulation age
nbyage_test_all$age[as.numeric(nbyage_test_all$age) > 15] <- 15

# Plot yearly nbyage
nbyage_test_all$age <- as.numeric(nbyage_test_all$age)
nbyage_test_all$model <- factor(nbyage_test_all$model, 
                                levels = c("high (1988-1999)", "all years", "low (2005-2019)"))

test_nbyage_plot <- ggplot(nbyage_test_all, aes(x=year, y=age)) +
  geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
  scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_y_continuous(breaks = seq(1, 15, 2), labels = c(seq(1, 13, 2), "15+")) +
  scale_x_discrete(breaks = seq(1988, 2022, 3)) +
  xlab(" ") + ylab("Age") +
  theme(legend.position = "none") +
  facet_wrap(~model, ncol = 1)
test_nbyage_plot

ggsave(filename = "plots/CEATTLE/cannibalism/Testing/dirichlet_nbyage.png", test_nbyage_plot,
       width=250, height=150, units="mm", dpi=300)


### Plot mortality ------------------------------------------------------------
# TODO: fix this - my own mortality plot code from ceattle_cannibalism.R is below - create faceted plot
# M2 <- extract_byage(run$quantities$M2, "multispecies", "M2")
# total_mortality <- M2 %>%
#   mutate(M1_M2 = M2 + rep(M1, length(years)))
# total_mortality$age <- as.integer(total_mortality$age)
# total_mortality$year <- as.integer(as.character(total_mortality$year))
# 
# mortality_plot <- ggplot(total_mortality, aes(y = age, x = year, zmin = 0, zmax = 1.5)) + 
#   geom_tile(aes(fill = M1_M2)) +
#   scale_y_continuous(expand = c(0, 0), breaks=c(1, 3, 5, 7, 9, 11, 13, 15)) + 
#   scale_x_continuous(expand = c(0, 0), breaks=c(1990, 1995, 2000, 2005, 2010, 2015, 2020)) + 
#   scale_fill_viridis(name = "M1 + M2", limits = c(0, 2.0), breaks = c(0.21, 2.0)) +
#   geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
#   coord_equal() +
#   ylab("Age") + xlab("Year") +
#   theme(panel.border = element_rect(colour = NA, fill = NA))
# 
# # Min, max, mean natural mortality by age, for ages 1-5
# M_byage <- total_mortality %>%
#   filter(age < 6) %>%
#   group_by(age) %>%
#   summarize(min = min(M1_M2), max = max(M1_M2), mean = mean(M1_M2))
# 
# ggsave(filename = "plots/CEATTLE/cannibalism/Testing/dirichlet_M.png", 
#        m_dirichlet, width=180, height = 180, units = "mm", dpi=300)
