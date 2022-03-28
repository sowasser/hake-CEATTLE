# Run CEATTLE for hake data with no diet data 
devtools::install_github("grantdadams/Rceattle@dev")
# library(Rceattle)

library(r4ss)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

# Run CEATTLE -----------------------------------------------------------------
# mydata <- Rceattle::read_data( file = "data/hake_from_ss.xlsx")
hake_intrasp <- Rceattle::read_data( file = "data/hake_intrasp_280310.xlsx")

intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                                 inits = NULL, # Initial parameters = 0
                                 file = NULL, # Don't save
                                 # debug = 1, # 1 = estimate, 0 = don't estimate
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 0, # Single species mode
                                 phase = "default")

# Check what all comes out of CEATTLE
ceattle_stuff <- intrasp_run$quantities

# Plot CEATTLE outputs
plot_biomass(Rceattle =  intrasp_run, incl_proj = FALSE)
plot_recruitment(Rceattle =  intrasp_run, add_ci = TRUE)
plot_catch(Rceattle =  intrasp_run, incl_proj = FALSE)
