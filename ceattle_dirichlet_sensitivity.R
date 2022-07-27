# Run CEATTLE for Pacific hake while testing different proportions of
# intraspecies predation, corrected with a Dirichlet multinomial distribution.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_220726.xlsx")

# # Run CEATTLE with the values as they are in the data file
# intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                  inits = NULL, # Initial parameters = 0
#                                  file = NULL, # Don't save
#                                  # debug = 1, # 1 = estimate, 0 = don't estimate
#                                  random_rec = FALSE, # No random recruitment
#                                  msmMode = 1, # Multispecies mode
#                                  phase = "default")


# Run CEATTLE with differing diet weight proportions --------------------------
# Read in different Dirichlet-corrected datasets
dirichlet_all <- read.csv("data/diet/Dirichlet/Dirichlet_all_years.csv")
dirichlet_90s <- read.csv("data/diet/Dirichlet/Dirichlet_90s.csv")
dirichlet_recent <- read.csv("data/diet/Dirichlet/Dirichlet_recent.csv")
intrasp_yearly <- read.csv("data/diet/diet_for_CEATTLE_yearly.csv")



# Adapt weight proportions to replace those in the excel file & run CEATTLE
run_ceattle <- function(df) {
  hake_intrasp$UobsWtAge <- df
  ceattle <- Rceattle::fit_mod(data_list = hake_intrasp,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # Multi-species mode
                               phase = "default")
  return(ceattle)
}

# Run low-cannibalism models
run_all <- run_ceattle(dirichlet_all)
run_90s <- run_ceattle(dirichlet_90s)
run_recent <- run_ceattle(dirichlet_recent)
run_tv <- run_ceattle(intrasp_yearly)

# Plot biomass in comparison to no diet diet run ---------------------------
years <- 1980:2022

# Pull out SSB & overall biomass from CEATTLE runs
ceattle_biomass <- function(run, name) {
  ssb <- (c(run$quantities$biomassSSB) * 2)
  biom <- c(run$quantities$biomass)
  wide <- as.data.frame(cbind(years, ssb, biom))
  colnames(wide) <- c("year", "SSB", "Total Biomass")
  all_biom <- melt(wide, id.vars = "year")
  colnames(all_biom)[2:3] <- c("type", name)
  
  return(all_biom)
}

test_biom <- cbind(ceattle_biomass(run_all, "Dirichlet - all years"),
                   ceattle_biomass(run_90s, "Dirichlet - 1991-1999"),
                   ceattle_biomass(run_recent, "Dirichlet - 2005-2019"),
                   ceattle_biomass(run_tv, "time-varying diet"))
test_biom <- test_biom[, c(1:3, 6, 9, 12)]

# Combine with run of CEATTLE using original data and plot
nodiet_biom <- read.csv("data/ceattle_nodiet_biom.csv")
colnames(nodiet_biom)[3] <- "no diet"

plot_biom <- function(df) {
  wide <- cbind(df, nodiet_biom[, 3])
  colnames(wide)[(ncol(wide))] <- c("no diet")
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
test_biom_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/test_dirichlet_biomass.png", test_biom_plot, 
       bg = "transparent", width=170, height=120, units="mm", dpi=300)


# Plot recruitment ------------------------------------------------------------
R_test_all <- cbind(c(run_all$quantities$R), c(run_90s$quantities$R), 
                    c(run_recent$quantities$R), c(run_tv$quantities$R))

nodiet_R <- read.csv("data/ceattle_nodiet_R.csv")
R_test_wide <- as.data.frame(cbind(years, R_test_all, nodiet_R))
colnames(R_test_wide) <- c("year",
                           "all years",
                           "1991-1999",
                           "2005-2019",
                           "time-varying diet",
                           "no diet")
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

ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/test_dirichlet_R.png", test_R_plot, 
       bg = "transparent", width=200, height=100, units="mm", dpi=300)


# # Calculate and plot difference btw no diet & each cannibalism run ------------
# nodiet_biom4 <- cbind(nodiet_biom[, 3], nodiet_biom[, 3],
#                       nodiet_biom[, 3], nodiet_biom[, 3])
# delta_biom_wide <- low_biom[44:86, c(3:6)] - nodiet_biom4
# delta_biom_wide <- cbind(years, delta_biom_wide)
# delta_biom <- melt(delta_biom_wide, id.vars = "years")
# 
# biom_difference <- ggplot(delta_biom2, aes(x=years, y=value)) +
#   geom_line(aes(color=variable)) +
#   scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.336) +  # colors to match biomass plot  
#   theme_sleek() +
#   ylab("Biomass (mt)") + xlab("year")
#   labs(color = "model") 
# biom_difference
# 
# ggsave(filename="plots/CEATTLE/nodietecies predation/Testing/low_nodiet_biom_difference.png", 
#        biom_difference, width=200, height=100, units="mm", dpi=300)


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

nodiet_nbyage <- read.csv("data/ceattle_nodiet_nbyage.csv")[, -4]
nodiet_nbyage <- cbind(nodiet_nbyage, rep("no diet"))
colnames(nodiet_nbyage)[4] <- "model"

nbyage_test_all <- rbind(extract_nbyage(run_all, "all years"),
                         extract_nbyage(run_90s, "1991-1999"),
                         extract_nbyage(run_recent, "2005-2019"),
                         extract_nbyage(run_tv, "time-varying diet"),
                         nodiet_nbyage)

# Set 15 as accumulation age
nbyage_test_all$age[as.numeric(nbyage_test_all$age) > 15] <- 15

# Calculate mean numbers at age & plot
nbyage_test_mean <- nbyage_test_all %>% group_by(age, model) %>%
  summarize(mean_number = mean(numbers))

# Change order of models to match other outputs
nbyage_test_mean$model <- factor(nbyage_test_mean$model, 
                                 levels = c("all years", 
                                            "1991-1999",
                                            "2005-2019",
                                            "time-varying diet",
                                            "no diet"))

test_nbyage_plot <- ggplot(nbyage_test_mean, aes(x=age, y=mean_number, fill=model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(labels = c(1:14, "15+")) +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  xlab("age") + ylab("numbers") 
test_nbyage_plot

ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_dirichlet_nbyage.png", test_nbyage_plot, 
       bg = "transparent", width=200, height=100, units="mm", dpi=300)


### Plot mortality ------------------------------------------------------------
m_all <- plot_mortality(Rceattle = run_all, type = 0, zlim = c(0, 3.2), title = "all years", maxage = 15) 
m_90s <- plot_mortality(Rceattle = run_90s, type = 0, zlim = c(0, 3.2), title = "1991-1999", maxage = 15) 
m_recent <- plot_mortality(Rceattle = run_recent, type = 0, zlim = c(0, 3.2), title = "2005-2019", maxage = 15) 
m_tv <- plot_mortality(Rceattle = run_tv, type = 0, zlim = c(0, 3.2), title = "time-varying diet", maxage = 15)


ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/M_dirichlet_all.png", m_all, dpi=300)
ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/M_dirichlet_90s.png", m_90s, dpi=300)
ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/M_dirichlet_recent.png", m_recent, dpi=300)
ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/M_dirichlet_recent.png", m_tv, dpi=300)
