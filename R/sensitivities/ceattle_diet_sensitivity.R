# Run CEATTLE for Pacific hake while testing different proportions of 
# cannibalism in their diet.

# devtools::install_github("grantdadams/Rceattle@dev")
library(reshape2)
library(dplyr)
library(purrr)
library(scales)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Read in prior M1 models
load("models/ms_priorM1.Rdata")
load("models/sensitivity/diet/run_wt05_prior.Rdata")
load("models/sensitivity/diet/run_wt10_prior.Rdata")
load("models/sensitivity/diet/run_wt50_prior.Rdata")
load("models/sensitivity/diet/run_wt75_prior.Rdata")

# Rceattle::plot_comp(run_wt05$model)
# Rceattle::plot_comp(run_wt10$model)
# Rceattle::plot_comp(run_wt50$model)
# Rceattle::plot_comp(run_wt75$model)

# Check fit of CEATTLE model --------------------------------------------------
sensitivity_fit <- rbind(cbind(model = "wt05", run_wt05_prior$fit),
                         cbind(model = "wt10", run_wt10_prior$fit),
                         cbind(model = "wt50", run_wt50_prior$fit),
                         cbind(model = "wt75", run_wt75_prior$fit))

sensitivity_summary <- list(run_wt05_prior$summary, run_wt10_prior$summary, 
                            run_wt50_prior$summary, run_wt75_prior$summary) %>%
  reduce(full_join, by = "component")
sensitivity_summary$NLL.x <- round(sensitivity_summary$NLL.x, digits = 1)
sensitivity_summary$NLL.y <- round(sensitivity_summary$NLL.y, digits = 1)
sensitivity_summary$NLL.x.x <- round(sensitivity_summary$NLL.x.x, digits = 1)
sensitivity_summary$NLL.y.y <- round(sensitivity_summary$NLL.y.y, digits = 1)
colnames(sensitivity_summary) <- c("component", "wt05", "wt10", "wt50", "wt75")


# Plot biomass & recruitment in comparison to original diet run ---------------
years <- 1980:2100
max_age <- max_age <- ms_priorM1$model$data_list$nages
models <- list(ms_priorM1$model, run_wt05_prior$model, run_wt10_prior$model, 
               run_wt50_prior$model, run_wt75_prior$model)
names <- c("Base Cannibalism", "Max 0.05", "Max 0.1", 
           "Max 0.5", "Max 0.75")

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

# Numbers-at-age
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

nbyage_test_all <- rbind(extract_byage(models[[1]]$quantities$NByage, names[1], "numbers"),
                         extract_byage(models[[2]]$quantities$NByage, names[2], "numbers"),
                         extract_byage(models[[3]]$quantities$NByage, names[3], "numbers"),
                         extract_byage(models[[4]]$quantities$NByage, names[4], "numbers"),
                         extract_byage(models[[5]]$quantities$NByage, names[5], "numbers"))

# Plot yearly nbyage
nbyage_test_all$age <- as.numeric(nbyage_test_all$age)
# nbyage_test_all$year <- as.numeric(nbyage_test_all$year)

test_nbyage_plot <- ggplot(nbyage_test_all, aes(x=year, y=age)) +
  geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
  scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_x_discrete(limits = c(as.character(1980:2022)), breaks = seq(1980, 2022, 3)) +
  xlab(" ") + ylab("Age") +
  theme(legend.position = "none") +
  facet_wrap(~model, ncol = 2, scales = "free_x")
test_nbyage_plot

# Calculate annual mean age
mean_nbyage <- nbyage_test_all %>%
  group_by(year, model) %>%
  summarize(value = weighted.mean(age, numbers)) %>%
  mutate(error = 0) %>%
  mutate(variable = "Mean Age")
mean_nbyage <- mean_nbyage[, c("year", "variable", "value", "error", "model")]

### Bring it all together! ----------------------------------------------------
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
  
  # Mean difference between models
  rel_change <- function(model1, model2, stat) {
    df1 <- all_popdy %>% 
      filter(model == model1 & variable == stat)
    df2 <- all_popdy %>%
      filter(model == model2) %>% filter(variable == stat) 
    mean_out <- mean((df1$value) - (df2$value))
    SEM <- sd((df1$value) - (df2$value)) / sqrt(length(range))
    percent <- mean(((df1$value - df2$value) / df1$value) * 100) 
    label <- paste(model1, "-", model2, ", ", stat)
    return(c(label, mean_out, SEM, percent))
  }
  
  rechange_all <- rbind(rel_change(names[2], names[1], "Spawning Biomass (Mt)"),
                        rel_change(names[2], names[1], "Total Biomass (Mt)"),
                        rel_change(names[2], names[1], "Recruitment (millions)"),
                        rel_change(names[3], names[1], "Spawning Biomass (Mt)"),
                        rel_change(names[3], names[1], "Total Biomass (Mt)"),
                        rel_change(names[3], names[1], "Recruitment (millions)"),
                        rel_change(names[4], names[1], "Spawning Biomass (Mt)"),
                        rel_change(names[4], names[1], "Total Biomass (Mt)"),
                        rel_change(names[4], names[1], "Recruitment (millions)"),
                        rel_change(names[5], names[1], "Spawning Biomass (Mt)"),
                        rel_change(names[5], names[1], "Total Biomass (Mt)"),
                        rel_change(names[5], names[1], "Recruitment (millions)"))
  
  all_popdy$model <- factor(all_popdy$model, 
                            levels = c(names[5], names[4], names[3], 
                                       names[2], names[1]))
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
    geom_line(aes(linetype = model)) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed")) +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) + 
    ylim(0, NA) + 
    xlim(1980, 2022) +
    ylab(" ") + xlab("Year") +
    labs(color = "Model", fill = "Model", linetype = "Model") +
    facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside") 
  
  return(list(all = all_popdy, change = rechange_all, plot = popdy_plot))
}

popdy_test <- plot_popdy(test_biom, R_test)

relative_change_test <- popdy_test$change

diet_popdy <- popdy_test$plot
diet_popdy

write.csv(popdy_test$all, file = "sensitivities/diet_popdy.csv", row.names = FALSE)

# ### Compare survey biomass estimate from CEATTLE to true values ---------------
# survey <- read.csv("data/assessment/survey_data.csv")
# survey <- cbind(survey, model = rep("SS3", length(survey$year)))
# 
# extract_srv <- function(run, name){
#   df <- data.frame(year = 1995:2019,
#                    biomass = run$quantities$srv_bio_hat,
#                    log_sd = run$quantities$srv_log_sd_hat,
#                    model = rep(name, length(1995:2019)))
#   return(df)
# }
# 
# srv_test <- rbind(extract_srv(run_wt05, "CEATTLE - 0.5% cannibalism"),
#                   extract_srv(run_wt10, "CEATTLE - 10% cannibalism"),
#                   extract_srv(run_wt50, "CEATTLE - 50% cannibalism"),
#                   extract_srv(run_wt80, "CEATTLE - 80% cannibalism"),
#                   survey)
# 
# test_survey_plot <- ggplot(srv_test, aes(x=year, y=biomass, color=model)) +
#   geom_line(alpha = 0.3) +
#   geom_point() +
#   # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
#   scale_color_viridis(discrete = TRUE, direction = -1) +
#   scale_fill_viridis(discrete = TRUE, direction = -1) +
#   xlab("year") + ylab("survey biomass") 
# test_survey_plot
# 
# ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_intrasp_survey.png", test_survey_plot, 
#        bg = "transparent", width=200, height=120, units="mm", dpi=300)


### Look at mortality-at-age timeseries ---------------------------------------
# Custom mortality function that outputs ggplot object ------------------------
mortality <- function(run, name) {
  M1 <- run$quantities$M1[1, 1, 1:max_age]
  M2 <- extract_byage(run$quantities$M2, name, "M2")
  total_mortality <- M2 %>%
    mutate(M1_M2 = M2 + rep(M1, length(years)))
  total_mortality$age <- as.integer(total_mortality$age)
  total_mortality$year <- as.integer(as.character(total_mortality$year))
  return(total_mortality)
}

M_all <- rbind(mortality(run_wt05_prior$model, names[2]),
               mortality(run_wt10_prior$model, names[3]),
               mortality(run_wt50_prior$model, names[4]),
               mortality(run_wt75_prior$model, names[5]))

max(M_all$M1_M2)  # check max M for plotting
M_test <- ggplot(M_all, aes(y = age, x = year, zmin = 0)) +
  geom_tile(aes(fill = M1_M2)) +
  scale_y_continuous(expand = c(0, 0), breaks=c(1, 5, 10, 15, 20)) +
  scale_x_continuous(expand = c(0, 0), limits = c(1980, 2022)) +
  scale_fill_viridis(name = "M1 + M2", limits = c(0, 2), breaks = c(0.21, 1, 2)) +
  geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  coord_equal() +
  ylab("Age") + xlab("Year") +
  facet_wrap(~model, ncol = 1)
M_test

# Min, max, mean natural mortality by age, for ages 1-5
M_byage <- M_all %>%
  filter(age < 6) %>%
  group_by(age, model) %>%
  summarize(min = min(M1_M2), max = max(M1_M2), mean = mean(M1_M2))

M1_all <- rbind(data.frame(model = "0.5%", 
                           mean = mean(run_wt05_prior$model$quantities$M1[1, 1, 1:15])),
                data.frame(model = "10%", 
                           mean = mean(run_wt10_prior$model$quantities$M1[1, 1, 1:15])),
                data.frame(model = "50%", 
                           mean = mean(run_wt50_prior$model$quantities$M1[1, 1, 1:15])),
                data.frame(model = "75%", 
                           mean = mean(run_wt75_prior$model$quantities$M1[1, 1, 1:15])))


### Save plots (when not experimenting) ---------------------------------------
# ggsave(filename = "plots/CEATTLE/cannibalism/Testing/sensitivity_nbyage.png", test_nbyage_plot,
#        width=220, height=210, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/Testing/sensitivity_meanage_popdy.png", popdy_test$plot,
#        width=140, height=170, units="mm", dpi=300)
# ggsave(filename = "plots/CEATTLE/cannibalism/Testing/sensitivity_M.png", M_test,
#        width=160, height = 250, units = "mm", dpi=300)
