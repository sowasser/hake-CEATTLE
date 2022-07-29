# Run CEATTLE with intraspecies-predation proportions calculated from diet 
# database going back to 1980.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

# Read in CEATTLE data from the excel file
hake_forecast <- Rceattle::read_data(file = "data/hake_forecast_220729.xlsx")

proj_run <- Rceattle::fit_mod(data_list = hake_forecast,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              # debug = 1, # 1 = estimate, 0 = don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 1, # Multi-species mode
                              phase = "default")

# Rceattle diagmostics
plot_index(proj_run)
plot_catch(proj_run)
plot_selectivity(proj_run)
plot_mortality(proj_run)
plot_indexresidual(proj_run)
plot_logindex(proj_run)
plot_recruitment(proj_run, add_ci = TRUE)
plot_comp(proj_run)
