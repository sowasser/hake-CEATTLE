# Run CEATTLE for hake data with no diet data 
devtools::install_github("grantdadams/Rceattle")
library(Rceattle)

library(r4ss)

# Run CEATTLE -----------------------------------------------------------------
mydata <- Rceattle::read_data( file = "data/hake_from_ss.xlsx")
grantsdata <- Rceattle::read_data( file = "data/From Grant/PacificHake/Data/2019PacificHake.xlsx")

# mydata$est_M1 <- c(0,0,0)
# mydata$estDynamics = 1
# mydata$fleet_control$proj_F_prop <-rep(0,2)  # changed from 7 in tutorial

ss_run <- Rceattle::fit_mod(data_list = grantsdata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            # debug = 1, # 1 = estimate, 0 = don't estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default")

plot_biomass(Rceattle =  ss_run)
plot_recruitment(Rceattle =  ss_run, add_ci = TRUE)
plot_catch(Rceattle =  ss_run, incl_proj = T)

# Compare estimates to assessment ---------------------------------------------
# Biomass from CEATTLE
biomass <- ss_run$quantities$biomass
biomass <- cbind(c(1966:2050), c(biomass))

# Diagnostics from r4ss
mydir <- file.path(file.path("hake_assessment/2020_Hake_Assessment"))

# read the model output and print diagnostic messages 
replist <- SS_output(dir = mydir, 
                     verbose = TRUE,
                     printstats = TRUE)

# plots the results
SS_plots(replist)

