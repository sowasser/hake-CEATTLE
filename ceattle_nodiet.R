# Run CEATTLE for hake data with no diet data 
# devtools::install_github("grantdadams/Rceattle")
library(Rceattle)

library(r4ss)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

# Run CEATTLE -----------------------------------------------------------------
# mydata <- Rceattle::read_data( file = "data/hake_from_ss.xlsx")
hake_nodiet <- Rceattle::read_data( file = "data/hake220203.xlsx")

# mydata$projyr <- 2022  # change to 2022 to match assessment
hake_nodiet$projyr <- 2022  # change to 2022 to match assessment

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

plot_biomass(Rceattle =  nodiet_run)
plot_recruitment(Rceattle =  nodiet_run, add_ci = TRUE)
plot_catch(Rceattle =  nodiet_run, incl_proj = F)

plot_ssb(Rceattle = nodiet_run)
plot_selectivity(Rceattle = nodiet_run)
plot_logindex(Rceattle = nodiet_run)

# Check what all comes out of CEATTLE
ceattle_stuff <- nodiet_run$quantities


# Compare spawning stock biomass & total biomass between SS3 and CEATTLE ------
# Pull out SSB & total biomass from CEATTLE & combine
ceattle_ssb <- c(nodiet_run$quantities$biomassSSB)
ceattle_biomass <- c(nodiet_run$quantities$biomass)

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

# Error from stock synthesis SSB
error <- c(rep(0, 114), ss_ssb_werror[, 3], rep(0, 57))

# Combine all together and plot
biom_wide <- cbind(ceattle_biom, ss_biom[, 3])
colnames(biom_wide)[4] <- "Stock Synthesis"
biom_noerror <- melt(biom_wide, id.vars = c("year", "type"))
biom <- cbind(biom_noerror, error)

biom_plot <- ggplot(biom, aes(x=year, y=value)) +
  geom_line(aes(color=variable, linetype=type)) +
  scale_linetype_manual(values=c("dashed", "solid")) +  # specify line types
  geom_errorbar(aes(ymin=value-error, ymax=value+error, color=variable), width=0, alpha=0.3) +
  scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +  # specify colors
  theme_sleek() +
  ylab("Biomass (mt)") +
  labs(color = "model")
biom_plot

ggsave(filename="plots/CEATTLE/allbiom_ss3_ceattle.png", biom_plot,
       width=200, height=100, units="mm", dpi=300)


# Compare recruitment between SS3 and CEATTLE ---------------------------------
ceattle_R <- c(nodiet_run$quantities$R)
ss_R <- read.table("data/assessment/recruitment.txt")

recruitment_wide <- as.data.frame(cbind(1966:2022, ceattle_R, ss_R[, 2]))
colnames(recruitment_wide) <- c("year", "CEATTLE", "Stock Synthesis")
recruitment <- melt(recruitment_wide, id.vars = "year")

# Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS is 0)
ss_1 <- cbind(1967:2022, rep("SS + 1", 56), ss_R[1:56, 2])
colnames(ss_1) <- c("year", "variable", "value")
recruitment <- rbind(recruitment, ss_1)
recruitment$value <- as.numeric(recruitment$value)
recruitment$year <- as.numeric(recruitment$year)

recruit_plot <- ggplot(recruitment, aes(x=year, y=value)) +
  geom_line(aes(color=variable, linetype=variable)) +
  scale_color_viridis(discrete = TRUE, begin = 0.15, end = 0.85) +  # specify colors
  scale_linetype_manual(values=c("solid", "dotted", "solid")) +
  theme_sleek() +
  ylab("Recruitment") +
  labs(color = "model")
recruit_plot

ggsave(filename="plots/CEATTLE/recruitment_ss3_ceattle.png", recruit_plot,
       width=200, height=100, units="mm", dpi=300)


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
