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


# Compare spawning stock biomass & total biomass between SS3 and CEATTLE ------
# Pull out SSB & total biomass from CEATTLE & combine
ceattle_ssb <- (c(nodiet_run$quantities$biomassSSB) * 2)
ceattle_biomass <- c(nodiet_run$quantities$biomass)

ceattle_biom_wide <- as.data.frame(cbind(years, ceattle_ssb, ceattle_biomass))
colnames(ceattle_biom_wide) <- c("year", "SSB", "Total Biomass")
ceattle_biom <- melt(ceattle_biom_wide, id.vars = "year")
colnames(ceattle_biom)[2:3] <- c("type", "CEATTLE no diet")

write.csv(ceattle_biom, "data/ceattle_nodiet_biom.csv", row.names = FALSE)


# Pull out SSB & total biomass from stock synthesis & combine, remove pre-1980
ss_ssb_werror <- read.table("data/assessment/ssb.txt")
ss_ssb <- ss_ssb_werror[15:57, 2]

ss_biomass <- read.table("data/assessment/biomass.txt")
ss_biom <- ss_biomass[15:57, 2]

ss_biom_wide <- as.data.frame(cbind(years, ss_ssb, ss_biom))
colnames(ss_biom_wide) <- c("year", "SSB", "total biomass")
ss_biom_all <- melt(ss_biom_wide, id.vars = "year")

# Combine all together and plot
biom_wide <- cbind(ceattle_biom, ss_biom_all[, 3])
colnames(biom_wide)[4] <- "Stock Synthesis"
biom <- melt(biom_wide, id.vars = c("year", "type"))

biom_plot <- ggplot(biom, aes(x=year, y=value)) +
  geom_line(aes(color=variable)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +  # specify colors
  theme_sleek() +
  ylab("Biomass (mt)") +
  labs(color = "model") +
  facet_wrap(~type, ncol = 1)
biom_plot

ggsave(filename="plots/CEATTLE/allbiom_ss3_ceattle.png", biom_plot,
       width=200, height=150, units="mm", dpi=300)


# Compare recruitment between SS3 and CEATTLE ---------------------------------
ceattle_R <- c(nodiet_run$quantities$R)
write.csv(ceattle_R, "data/ceattle_nodiet_R.csv", row.names = FALSE)

ss_R <- read.table("data/assessment/recruitment.txt")[15:57,]

recruitment_wide <- as.data.frame(cbind(years, ceattle_R, ss_R[, 2]))
colnames(recruitment_wide) <- c("year", "CEATTLE", "Stock Synthesis")
recruitment <- melt(recruitment_wide, id.vars = "year")

# Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS is 0)
ss_1 <- cbind(1981:2022, rep("SS + 1", (length(years)-1)), ss_R[1:42, 2])
colnames(ss_1) <- c("year", "variable", "value")
recruitment <- rbind(recruitment, ss_1)
recruitment$value <- as.numeric(recruitment$value)
recruitment$year <- as.numeric(recruitment$year)

recruit_plot <- ggplot(recruitment, aes(x=year, y=value)) +
  geom_line(aes(color=variable, linetype=variable)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +  # specify colors
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
