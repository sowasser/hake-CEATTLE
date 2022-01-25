# Run CEATTLE for hake data with no diet data 
library(devtools)
# library(Rceattle)
install_github("grantdadams/Rceattle")

mydata <- Rceattle::read_data( file = "data/input_nodiet.xlsx")
mydata$est_M1 <- c(0,0,0)
mydata$estDynamics = 0

mydata$fleet_control$proj_F_prop <-rep(0,3)  # changed from 7 in tutorial
ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            debug = FALSE, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default")




# Grant's suggestions
################################################
# Data
################################################
# Read the data in
mydata_pollock_fixed <- Rceattle::read_data( file = "GOA_18.5.1_pollock_single_species_1970-2018.xlsx")
mydata_pollock <-  mydata_pollock_fixed
mydata_pollock$estDynamics = 0

#######################################
# Mod 1 - Fix n-at-age
#######################################
mydata_pollock_fixed$estDynamics = 1

pollock_fixed <- Rceattle::fit_mod(data_list = mydata_pollock_fixed,
                                   inits = NULL, # Initial parameters = 0
                                   file = NULL, # Don't save
                                   debug = 1, # Don't estimate (no need to)
                                   msmMode = 0, # Single species mode
                                   silent = TRUE,
                                   phase = "default")

#######################################
# Mod 2 - Base
#######################################

pollock_base <- Rceattle::fit_mod(data_list = mydata_pollock,
                                  inits = NULL, # Initial parameters = 0
                                  file = NULL, # Don't save
                                  debug = 0, # Estimate
                                  msmMode = 0, # Single species mode
                                  silent = TRUE,
                                  phase = "default")

# Compare

plot_biomass(list(pollock_fixed, pollock_base), model_names = c("Fixed","Est"))