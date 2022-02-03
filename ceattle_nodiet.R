# Run CEATTLE for hake data with no diet data 
devtools::install_github("grantdadams/Rceattle")
library(Rceattle)

library(r4ss)
library(reshape2)
library(ggplot2)
library(ggsidekick)

# Run CEATTLE -----------------------------------------------------------------
mydata <- Rceattle::read_data( file = "data/hake_from_ss.xlsx")
GAdata <- Rceattle::read_data( file = "data/From Grant/PacificHake/Data/2019PacificHake.xlsx")

GAdata$projyr <- 2022  # change to 2022 to match assessment
mydata$projyr <- 2022  # change to 2022 to match assessment

ss_run_GA <- Rceattle::fit_mod(data_list = GAdata,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                             # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 0, # Single species mode
                               phase = "default")

ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            # debug = 1, # 1 = estimate, 0 = don't estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default")

plot_biomass(Rceattle =  ss_run)
plot_recruitment(Rceattle =  ss_run, add_ci = TRUE)
plot_catch(Rceattle =  ss_run, incl_proj = T)


# Compare virgin spawning stock biomass between SS3 and CEATTLE ---------------
# Biomass from CEATTLE
ceattle_biomass <- c(ss_run$quantities$biomass)
ceattle_biomass <- cbind(1966:2022, ceattle_biomass)

# Pull out just biomass from stock synthesis run
ss_biomass <- read.table("data/assessment/ssb.txt")

# Combine and plot
biomass_wide <- as.data.frame(cbind(ceattle_biomass, ss_biomass[, 2]))
colnames(biomass_wide) <- c("year", "CEATTLE", "SS3")
biomass <- melt(biomass_wide, id.vars = "year")

ssb_plot <- ggplot(biomass, aes(x=year, y=value)) +
  geom_line(aes(color=variable), size=1) +
  # scale_color_viridis(discrete = TRUE) +  # invert colors
  theme_sleek() +
  ylab("biomass") +
  labs(color = "model")
ssb_plot


# Run r4ss diagnostics on stock synthesis model -------------------------------
# Diagnostics from r4ss
mydir <- file.path(file.path("hake_assessment/2020_Hake_Assessment"))

# read the model output and print diagnostic messages 
replist <- SS_output(dir = mydir, 
                     verbose = TRUE,
                     printstats = TRUE)

# plots the results
SS_plots(replist)
