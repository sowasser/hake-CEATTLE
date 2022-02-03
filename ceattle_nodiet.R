# Run CEATTLE for hake data with no diet data 
devtools::install_github("grantdadams/Rceattle")
library(Rceattle)

library(r4ss)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

# Run CEATTLE -----------------------------------------------------------------
# mydata <- Rceattle::read_data( file = "data/hake_from_ss.xlsx")
GAdata <- Rceattle::read_data( file = "data/From Grant/PacificHake/Data/2019PacificHake.xlsx")

# mydata$projyr <- 2022  # change to 2022 to match assessment
GAdata$projyr <- 2022  # change to 2022 to match assessment

# ss_run <- Rceattle::fit_mod(data_list = mydata,
#                             inits = NULL, # Initial parameters = 0
#                             file = NULL, # Don't save
#                             # debug = 1, # 1 = estimate, 0 = don't estimate
#                             random_rec = FALSE, # No random recruitment
#                             msmMode = 0, # Single species mode
#                             phase = "default")

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


# Compare spawning stock biomass & total biomass between SS3 and CEATTLE ------
# Pull out SSB & total biomass from CEATTLE & combine
ceattle_ssb <- c(ss_run_GA$quantities$biomassSSB)
ceattle_biomass <- c(ss_run_GA$quantities$biomass)

ceattle_biom_wide <- as.data.frame(cbind(1966:2022, ceattle_ssb, ceattle_biomass))
colnames(ceattle_biom_wide) <- c("year", "SSB", "Total Biomass")
ceattle_biom <- melt(ceattle_biom_wide, id.vars = "year")
colnames(ceattle_biom)[2:3] <- c("type", "CEATTLE")

# Pull out SSB & total biomass from stock synthesis & combine
ss_ssb_werror <- read.table("data/assessment/ssb.txt")
ss_ssb <- ss_ssb_werror[, 2]

ss_biomass <- read.table("data/assessment/biomass.txt", header = TRUE)
ss_biomass2 <- ss_biomass %>% filter(Beg.Mid == "B")
ss_biomass3 <- ss_biomass2[3:59, 13:33]
ss_biomass <- rowSums(ss_biomass3)

ss_biom_wide <- as.data.frame(cbind(1966:2022, ss_ssb, ss_biomass))
colnames(ss_biom_wide) <- c("year", "SSB", "total biomass")
ss_biom <- melt(ss_biom_wide, id.vars = "year")

# Combine all together and plot
biom_wide <- cbind(ceattle_biom, ss_biom[, 3])
colnames(biom_wide)[4] <- "Stock Synthesis"
biom <- melt(biom_wide, id.vars = c("year", "type"))

biom_plot <- ggplot(biom, aes(x=year, y=value)) +
  geom_line(aes(color=variable, linetype=type)) +
  scale_linetype_manual(values=c("dashed", "solid")) +  # specify line types
  scale_color_viridis(discrete = TRUE, begin = 0.25, end = 0.75) +  # specify colors
  theme_sleek() +
  ylab("Biomass") +
  labs(color = "model")
biom_plot

ggsave(filename="plots/CEATTLE/allbiom_ss3_ceattle.png", biom_plot,
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
