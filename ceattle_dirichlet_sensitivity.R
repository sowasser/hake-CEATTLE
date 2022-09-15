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


### Run CEATTLE with differing diet weight proportions ------------------------
# Read in different Dirichlet-corrected datasets
dirichlet_90s <- read.csv("data/diet/Dirichlet/Dirichlet_90s.csv")
dirichlet_recent <- read.csv("data/diet/Dirichlet/Dirichlet_recent.csv")


# Adapt weight proportions to replace those in the excel file & run CEATTLE
run_ceattle <- function(df, start, end, proj) {
  hake_intrasp$UobsWtAge <- df
  hake_intrasp$styr <- start
  hake_intrasp$endyr <- end
  hake_intrasp$projyr <- proj
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
run_all <- run_ceattle(hake_intrasp$UobsWtAge, 1980, 2019, 2022)
run_90s <- run_ceattle(dirichlet_90s, 1991, 1999, 1999)
run_recent <- run_ceattle(dirichlet_recent, 2005, 2019, 2019)

# Run CEATTLE with no diet
nodiet_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                                inits = NULL, # Initial parameters = 0
                                file = NULL, # Don't save
                                # debug = 1, # 1 = estimate, 0 = don't estimate
                                random_rec = FALSE, # No random recruitment
                                msmMode = 0, # Single-species mode
                                phase = "default")


### Plot population dynamics --------------------------------------------------
ceattle_popdy <- function(run, name, years) {
  ssb <- (c(run$quantities$biomassSSB) * 2)
  biom <- c(run$quantities$biomass)
  wide <- as.data.frame(cbind(years, ssb, biom))
  colnames(wide) <- c("year", "SSB", "Total Biomass")
  all_biom <- melt(wide, id.vars = "year")
  all_biom2 <- cbind(all_biom, model = rep(name, length(all_biom$year)))  
  colnames(all_biom2)[2:3] <- c("type", "value")
  
  recruitment <- cbind(year = years,
                       type = rep("Recruitment"),
                       value = c(run$quantities$R),
                       model = name)
  
  popdy <- rbind(all_biom2, recruitment)
  
  return(popdy)
}

popdy <- rbind(ceattle_popdy(run_all, "all years", 1980:2022), 
               ceattle_popdy(run_90s, "1991-1999", 1991:1999), 
               ceattle_popdy(run_recent, "2005-2019", 2005:2019),
               ceattle_popdy(nodiet_run, "single-species", 1980:2022))

popdy$year <- as.numeric(popdy$year)
popdy$value <- as.numeric(popdy$value)
popdy$model <- factor(popdy$model,
                      levels = c("all years", "single-species", "1991-1999", "2005-2019"))

popdy_plot <- ggplot(popdy, aes(x=year, y=value, color = model, fill = model)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
  ylab(" ") +
  labs(color = "model") +
  facet_wrap(~type, ncol = 1, scales = "free_y")
popdy_plot


### Plot mortality ------------------------------------------------------------
m_all <- plot_mortality(Rceattle = run_all, type = 0, title = "all years", maxage = 15) 
m_90s <- plot_mortality(Rceattle = run_90s, type = 0, title = "1991-1999", maxage = 15) 
m_recent <- plot_mortality(Rceattle = run_recent, type = 0, title = "2005-2019", maxage = 15) 

ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/M_dirichlet_all.png", m_all, dpi=300)
ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/M_dirichlet_90s.png", m_90s, dpi=300)
ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/M_dirichlet_recent.png", m_recent, dpi=300)
