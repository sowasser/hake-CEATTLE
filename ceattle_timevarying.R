# Run CEATTLE for Pacific hake with time-varying predation and compare with the
# base model.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_220628.xlsx")

# # Run CEATTLE with the values as they are in the data file
# intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                  inits = NULL, # Initial parameters = 0
#                                  file = NULL, # Don't save
#                                  # debug = 1, # 1 = estimate, 0 = don't estimate
#                                  random_rec = FALSE, # No random recruitment
#                                  msmMode = 1, # Multispecies mode
#                                  phase = "default")


# Run CEATTLE with time-varying diet ------------------------------------------
intrasp_yearly <- read.csv("data/diet/diet_for_CEATTLE_yearly.csv")

df <- hake_intrasp$UobsWtAge
hake_intrasp$UobsWtAge <- intrasp_yearly
ceattle <- Rceattle::fit_mod(data_list = hake_intrasp,
                             inits = NULL, # Initial parameters = 0
                             file = NULL, # Don't save
                             # debug = 1, # 1 = estimate, 0 = don't estimate
                             random_rec = FALSE, # No random recruitment
                             msmMode = 1, # Multi-species mode
                             phase = "default")

years <- 1980:2022

# Pull out SSB & overall biomass from CEATTLE runs ----------------------------
ceattle_biomass <- function(run, name) {
  ssb <- (c(run$quantities$biomassSSB) * 2)
  biom <- c(run$quantities$biomass)
  wide <- as.data.frame(cbind(years, ssb, biom))
  colnames(wide) <- c("year", "SSB", "Total Biomass")
  all_biom <- melt(wide, id.vars = "year")
  colnames(all_biom)[2:3] <- c("type", name)
  
  return(all_biom)
}

test_biom <- ceattle_biomass(ceattle, "time-varying")

# Combine with run of CEATTLE using original data and plot
intrasp_biom <- read.csv("data/ceattle_intrasp_biomass.csv")
colnames(intrasp_biom)[3] <- "CEATTLE - unweighted"

plot_biom <- function(df) {
  wide <- cbind(test_biom, intrasp_biom[, 3])
  colnames(wide)[(ncol(wide))] <- c("static")
  biom <- melt(wide, id.vars = c("year", "type"))
  
  plot <- ggplot(biom, aes(x=year, y=value)) +
    geom_line(aes(color=variable)) +
    scale_color_viridis(discrete = TRUE, direction = -1) +  
    ylab("Biomass (mt)") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 1)
  
  return(plot)
}

# Low-cannibalism plot
test_biom_plot <- plot_biom(test_biom)
test_biom_plot  # no difference!!!

# ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/test_dirichlet_biomass.png", test_biom_plot, 
#        bg = "transparent", width=170, height=120, units="mm", dpi=300)


# Plot recruitment ------------------------------------------------------------
R_test_all <- c(ceattle$quantities$R)

intrasp_R <- read.csv("data/ceattle_intrasp_R.csv")
R_test_wide <- as.data.frame(cbind(years, R_test_all, intrasp_R))
colnames(R_test_wide) <- c("year",
                           "time-varying",
                           "static")
R_test <- melt(R_test_wide, id.vars = "year")

plot_R <- function(df) {
  df$value <- as.numeric(df$value)
  df$year <- as.numeric(df$year)
  
  plot <- ggplot(df, aes(x=year, y=value)) +
    geom_line(aes(color=variable)) +
    scale_color_viridis(discrete = TRUE, direction = -1) + 
    ylab("Recruitment") +
    labs(color = "model")
  
  return(plot)
}

test_R_plot <- plot_R(R_test)
test_R_plot

# ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/test_dirichlet_R.png", test_R_plot, 
#        bg = "transparent", width=200, height=100, units="mm", dpi=300)


# Numbers-at-age for each model run -------------------------------------------
# Read in data from no diet CEATTLE run
extract_nbyage <- function(run, name) {
  df <- as.data.frame(as.table(run$quantities$NByage))
  
  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(1980:2022)
  colnames(df) <- c("age", "year", "numbers")
  
  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"
  
  return(df)
}

intrasp_nbyage <- read.csv("data/ceattle_intrasp_nbyage.csv")[, -4]
intrasp_nbyage <- cbind(intrasp_nbyage, rep("static"))
colnames(intrasp_nbyage)[4] <- "model"

nbyage_test_all <- rbind(extract_nbyage(ceattle, "time-varying"),
                         intrasp_nbyage)

# Set 15 as accumulation age
nbyage_test_all$age[as.numeric(nbyage_test_all$age) > 15] <- 15

# Calculate mean numbers at age & plot
nbyage_test_mean <- nbyage_test_all %>% group_by(age, model) %>%
  summarize(mean_number = mean(numbers))

# Change order of models to match other outputs
nbyage_test_mean$model <- factor(nbyage_test_mean$model, 
                                 levels = c("time-varying", "static")

test_nbyage_plot <- ggplot(nbyage_test_mean, aes(x=age, y=mean_number, fill=model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(labels = c(1:14, "15+")) +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  xlab("age") + ylab("numbers") 
test_nbyage_plot

# ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_dirichlet_nbyage.png", test_nbyage_plot, 
#        bg = "transparent", width=200, height=100, units="mm", dpi=300)
