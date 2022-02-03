# Run CEATTLE for hake data with no diet data 
devtools::install_github("grantdadams/Rceattle")
library(Rceattle)

library(r4ss)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsidekick)

# Run CEATTLE -----------------------------------------------------------------
mydata <- Rceattle::read_data( file = "data/hake_from_ss.xlsx")
GAdata <- Rceattle::read_data( file = "data/From Grant/PacificHake/Data/2019PacificHake.xlsx")

mydata$projyr <- 2022  # change to 2022 to match assessment
GAdata$projyr <- 2022  # change to 2022 to match assessment

ss_run <- Rceattle::fit_mod(data_list = mydata,
                            inits = NULL, # Initial parameters = 0
                            file = NULL, # Don't save
                            # debug = 1, # 1 = estimate, 0 = don't estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = "default")

ss_run_GA <- Rceattle::fit_mod(data_list = GAdata,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 0, # Single species mode
                               phase = "default")

plot_biomass(Rceattle =  ss_run_GA)
plot_recruitment(Rceattle =  ss_run_GA, add_ci = TRUE)
plot_catch(Rceattle =  ss_run_GA, incl_proj = F)

plot_ssb(Rceattle = ss_run_GA)
plot_selectivity(Rceattle = ss_run_GA)
plot_logindex(Rceattle = ss_run_GA)


# Compare spawning stock biomass between SS3 and CEATTLE ----------------------
ceattle_ssb <- c(ss_run_GA$quantities$biomassSSB)
ceattle_ssb <- cbind(1966:2022, ceattle_ssb)

# Pull out just biomass from stock synthesis run
ss_ssb <- read.table("data/assessment/ssb.txt")

# Combine and plot
ssb_wide <- as.data.frame(cbind(ceattle_ssb, ss_ssb[, 2]))
colnames(ssb_wide) <- c("year", "CEATTLE", "SS3")
ssb <- melt(ssb_wide, id.vars = "year")

ssb_plot <- ggplot(ssb, aes(x=year, y=value)) +
  geom_line(aes(color=variable), size=1) +
  # scale_color_viridis(discrete = TRUE) +  # invert colors
  theme_sleek() +
  ylab("spawning stock biomass") +
  labs(color = "model")
ssb_plot

ggsave(filename="plots/CEATTLE/SSB_ss3_ceattle.png", ssb_plot,
       width=200, height=100, units="mm", dpi=300)


# Compare biomass between SS3 and CEATTLE -------------------------------------
# Biomass from CEATTLE
ceattle_biomass <- c(ss_run_GA$quantities$biomass)
ceattle_biomass <- cbind(1966:2022, ceattle_biomass)

ss_biomass <- read.table("data/assessment/biomass.txt", header = TRUE)
ss_biomass2 <- ss_biomass %>% filter(Beg.Mid == "B")
ss_biomass3 <- ss_biomass2[3:59, 13:33]
ss_biom <- rowSums(ss_biomass3)

biomass_wide <- as.data.frame(cbind(ceattle_biomass, ss_biom))
colnames(biomass_wide) <- c("year", "CEATTLE", "SS3")
biomass <- melt(biomass_wide, id.vars = "year")

biomass_plot <- ggplot(biomass, aes(x=year, y=value)) +
  geom_line(aes(color=variable), size=1) +
  # scale_color_viridis(discrete = TRUE) +  # invert colors
  theme_sleek() +
  ylab("biomass") +
  labs(color = "model")
biomass_plot

ggsave(filename="plots/CEATTLE/biomass_ss3_ceattle.png", biomass_plot,
       width=200, height=100, units="mm", dpi=300)


# Run r4ss diagnostics on stock synthesis model -------------------------------
# Diagnostics from r4ss
mydir <- file.path(file.path("hake_assessment/2020_Hake_Assessment"))

# read the model output and print diagnostic messages 
replist <- SS_output(dir = mydir, 
                     verbose = TRUE,
                     printstats = TRUE)

# plots the results
SS_plots(replist)
