# Run CEATTLE for hake data with no diet data 
# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)

library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

hake_intrasp <- Rceattle::read_data( file = "data/hake_intrasp_220401.xlsx")

# # Run CEATTLE with the values as they are in the data file
# intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                  inits = NULL, # Initial parameters = 0
#                                  file = NULL, # Don't save
#                                  # debug = 1, # 1 = estimate, 0 = don't estimate
#                                  random_rec = FALSE, # No random recruitment
#                                  msmMode = 0, # Single species mode
#                                  phase = "default")


# Run CEATTLE with differing diet weight proportions --------------------------
# Set different diet weight proportion distributions
wt05 <- c(0.0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
wt10 <- c(0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
wt30 <- c(0.0, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)
wt50 <- c(0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
wt80 <- c(0.0, 0.16, 0.24, 0.32, 0.40, 0.48, 0.56, 0.64, 0.72, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)

zeros <- rep(0, 19)

# Adapt weight proportions to replace those in the excel file & run CEATTLE
run_ceattle <- function(wt, df) {
  new_wt <- c()
  for(i in wt) {
    new_wt <- append(new_wt, c(i, zeros))
  }
  df$UobsWtAge$Stomach_proportion_by_weight <- wt
  ceattle <- Rceattle::fit_mod(data_list = df,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 0, # Single species mode
                               phase = "default")
  return(ceattle)
}

wt05_run <- run_ceattle(wt05, hake_intrasp)
wt10_run <- run_ceattle(wt10, hake_intrasp)
wt30_run <- run_ceattle(wt30, hake_intrasp)
wt50_run <- run_ceattle(wt50, hake_intrasp)
wt80_run <- run_ceattle(wt80, hake_intrasp)


# Check what all comes out of CEATTLE
# ceattle_stuff <- intrasp_run$quantities

years <- 1980:2022


# Compare intraspecies CEATTLE runs to no diet & assessment -------------------
# Pull out SSB & overall biomass from CEATTLE runs
ceattle_biomass <- function(run, name) {
  ssb <- (c(run$quantities$biomassSSB) * 2)
  biom <- c(run$quantities$biomass)
  wide <- as.data.frame(cbind(years, ssb, biom))
  colnames(wide) <- c("year", "SSB", "Total Biomass")
  all_biom <- melt(wide, id.vars = "year")
  colnames(all_biom)[2:3] <- c("type", name)
  
  return(all_biom)
}

all_test <- cbind(ceattle_biomass(wt05_run, "CEATTLE - 0.5% cannibalism"), 
                  ceattle_biomass(wt10_run, "CEATTLE - 10% cannibalism"),
                  ceattle_biomass(wt30_run, "CEATTLE - 30% cannibalism"),
                  ceattle_biomass(wt50_run, "CEATTLE - 50% cannibalism"),
                  ceattle_biomass(wt80_run, "CEATTLE - 80% cannibalism"))
all_test <- all_test[, c(1:3, 6 ,9, 12, 15)]

# Read in no diet data
nodiet_biom <- read.csv("data/ceattle_nodiet_biom.csv")
colnames(nodiet_biom)[3] <- "CEATTLE - no diet"

# Pull out SSB & total biomass from stock synthesis & combine, remove pre-1980
ss_ssb_werror <- read.table("data/assessment/ssb.txt")
ss_ssb <- ss_ssb_werror[15:57, 2]

ss_biomass <- read.table("data/assessment/biomass.txt")
ss_biom <- ss_biomass[15:57, 2]

ss_biom_wide <- as.data.frame(cbind(years, ss_ssb, ss_biom))
colnames(ss_biom_wide) <- c("year", "SSB", "total biomass")
ss_biom_all <- melt(ss_biom_wide, id.vars = "year")

# Combine all together and plot
biom_wide <- cbind(all_test, nodiet_biom[, 3], ss_biom_all[, 3])
colnames(biom_wide)[8:9] <- c("CEATTLE - no diet", "Stock Synthesis")
biom <- melt(biom_wide, id.vars = c("year", "type"))

biom_plot <- ggplot(biom, aes(x=year, y=value)) +
  geom_line(aes(color=variable)) +
  scale_color_viridis(discrete = TRUE) +  
  theme_sleek() +
  ylab("Biomass (mt)") +
  labs(color = "model") +
  facet_wrap(~type, ncol = 1)
biom_plot

ggsave(filename="plots/CEATTLE/intrasp_test_biomass.png", biom_plot,
       width=200, height=150, units="mm", dpi=300)

# Compare recruitment
ceattle_R <- cbind(c(wt05_run$quantities$R), c(wt10_run$quantities$R), 
                   c(wt30_run$quantities$R), c(wt50_run$quantities$R),
                   c(wt80_run$quantities$R))

nodiet_R <- read.csv("data/ceattle_nodiet_R.csv")

ss_R <- read.table("data/assessment/recruitment.txt")[15:57,]

recruitment_wide <- as.data.frame(cbind(years, ceattle_R, nodiet_R))
colnames(recruitment_wide) <- c("year", 
                                "CEATTLE - 0.5% cannibalism", 
                                "CEATTLE - 10% cannibalism", 
                                "CEATTLE - 30% cannibalism", 
                                "CEATTLE - 50% cannibalism",
                                "CEATTLE - 80% cannibalism",
                                "CEATTLE - no diet")
recruitment <- melt(recruitment_wide, id.vars = "year")

# Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS is 0)
ss_1 <- cbind(1981:2022, rep("SS + 1", (length(years)-1)), ss_R[1:42, 2])
colnames(ss_1) <- c("year", "variable", "value")

recruitment <- rbind(recruitment, ss_1)
recruitment$value <- as.numeric(recruitment$value)
recruitment$year <- as.numeric(recruitment$year)
recruit_plot <- ggplot(recruitment, aes(x=year, y=value)) +
  geom_line(aes(color=variable)) +
  scale_color_viridis(discrete = TRUE) + 
  theme_sleek() +
  ylab("Recruitment") +
  labs(color = "model")
recruit_plot

ggsave(filename="plots/CEATTLE/intrasp_test_recruit.png", recruit_plot,
       width=200, height=100, units="mm", dpi=300)
