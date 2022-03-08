# Run CEATTLE for hake data with no diet data 
devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)

library(r4ss)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

# Run CEATTLE -----------------------------------------------------------------
# mydata <- Rceattle::read_data( file = "data/hake_from_ss.xlsx")
hake_nodiet <- Rceattle::read_data( file = "data/hake_singlesp_220307.xlsx")

# ss_run <- Rceattle::fit_mod(data_list = mydata,
#                             inits = NULL, # Initial parameters = 0
#                             file = NULL, # Don't save
#                             # debug = 1, # 1 = estimate, 0 = don't estimate
#                             random_rec = FALSE, # No random recruitment
#                             msmMode = 0, # Single species mode
#                             phase = "default")

nodiet_run <- Rceattle::fit_mod(data_list = hake_nodiet,
                                inits = NULL, # Initial parameters = 0
                                file = NULL, # Don't save
                                # debug = 1, # 1 = estimate, 0 = don't estimate
                                random_rec = FALSE, # No random recruitment
                                msmMode = 0, # Single species mode
                                phase = "default")

# Check what all comes out of CEATTLE
ceattle_stuff <- nodiet_run$quantities