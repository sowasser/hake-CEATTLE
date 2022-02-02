# Install necessary packages --------------------------------------------------
# install.packages("devtools")

# Install TMB and rtools https://cran.r-project.org/bin/windows/Rtools/
# Instructions here for non-pc: https://github.com/kaskr/adcomp/wiki/Download
# install.packages('TMB', type = 'source')
TMB::runExample(all = TRUE)  # see if TMB works

# Install Rceattle - https://github.com/grantdadams/Rceattle
devtools::install_github("grantdadams/Rceattle")
library(Rceattle)

# Single-species implementation -----------------------------------------------

################################################
# Data
################################################
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2030

# Write data to excel
Rceattle::write_data(data_list = BS2017SS, file = "data/BS2017SS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata <- Rceattle::read_data( file = "data/BS2017SS.xlsx")
mydata$est_M1 <- c(0,0,0)


################################################
# Estimation
################################################
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
mydata$fleet_control$proj_F_prop <-rep(0,7)
ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                          # debug = FALSE, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default")
                          # silent = TRUE) # now getting an "unused argument" error here

# Estimate pollock M
mydata_M <- mydata
mydata_M$est_M1 <- c(1,0,0)
ss_run_M <- Rceattle::fit_mod(data_list = mydata_M,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                            # debug = FALSE, # Estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = "default" )
                            # silent = TRUE )  # now getting an "unused argument" error here
                            # recompile = FALSE)  # get an "unused argument" error if this is included

# The you can plot the model results using using
plot_biomass(Rceattle =  list(ss_run, ss_run_M))
plot_recruitment(Rceattle =  ss_run, add_ci = TRUE)
plot_catch(Rceattle =  ss_run, incl_proj = T)
