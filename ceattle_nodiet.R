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
hake_nodiet <- Rceattle::read_data( file = "data/hake_singlesp_220401.xlsx")

nodiet_run <- Rceattle::fit_mod(data_list = hake_nodiet,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                            # debug = 1, # 1 = estimate, 0 = don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = "default")

# Check what all comes out of CEATTLE
ceattle_stuff <- nodiet_run$quantities

years <- 1980:2022

# Pull out SSB & total biomass from CEATTLE & combine -------------------------
nodiet_ssb <- (c(nodiet_run$quantities$biomassSSB) * 2)
nodiet_biomass <- c(nodiet_run$quantities$biomass)

nodiet_biom_wide <- as.data.frame(cbind(years, nodiet_ssb, nodiet_biomass))
colnames(nodiet_biom_wide) <- c("year", "SSB", "Total Biomass")
nodiet_biom <- melt(nodiet_biom_wide, id.vars = "year")
colnames(nodiet_biom)[2:3] <- c("type", "CEATTLE no diet")

write.csv(nodiet_biom, "data/ceattle_nodiet_biom.csv", row.names = FALSE)


# Compare recruitment between SS3 and CEATTLE ---------------------------------
nodiet_R <- c(nodiet_run$quantities$R)
write.csv(nodiet_R, "data/ceattle_nodiet_R.csv", row.names = FALSE)


# # Run r4ss diagnostics on stock synthesis model -------------------------------
# # Diagnostics from r4ss
# mydir <- file.path(file.path("hake_assessment/2020_Hake_Assessment"))
# 
# # read the model output and print diagnostic messages 
# replist <- SS_output(dir = mydir, 
#                      verbose = TRUE,
#                      printstats = TRUE)
# 
# # plots the results
# SS_plots(replist)
