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
#                                  msmMode = 1, # Multispecies mode
#                                  phase = "default")


# Run CEATTLE with differing diet weight proportions --------------------------
# Set different diet weight proportion distributions
wt001 <- c(0.0, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001)
wt005 <- c(0.0, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005)
wt01 <- c(0.0, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
wt05 <- c(0.0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)

# # Plot stomach contents curves
# # Pull out data from base intrasp run
# wts <- hake_intrasp$UobsWtAge$Stomach_proportion_by_weight
# wts_short <- wts[wts != 0]
# prop <- as.data.frame(cbind(1:15, wt001, wt005, wt01, wt05, wts_short))
# colnames(prop)[c(1, 6)] <- c("age", "Grant's")
# prop_all <- melt(prop, id.vars = "age")
# 
# stomach_props <- ggplot(prop_all, aes(x=age, y=value, fill=variable)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_sleek() +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   ylab("stomach proportion")
# stomach_props


# Adapt weight proportions to replace those in the excel file & run CEATTLE
run_ceattle <- function(wt, df) {
  new_wt <- c()
  for(i in wt) {
    new_wt <- append(new_wt, c(i, rep(0, 14)))
  }
  df$UobsWtAge$Stomach_proportion_by_weight <- new_wt
  ceattle <- Rceattle::fit_mod(data_list = df,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # Multi-species mode
                               phase = "default")
  return(ceattle)
}

# Run low-cannibalism models
run_wt001 <- run_ceattle(wt001, hake_intrasp)
run_wt005 <- run_ceattle(wt005, hake_intrasp)
run_wt01 <- run_ceattle(wt01, hake_intrasp)
run_wt05 <- run_ceattle(wt05, hake_intrasp)

# Check what all comes out of CEATTLE
# ceattle_stuff <- run_wt001$quantities


# Plot biomass in comparison to no diet & asssessment -------------------------
years <- 1980:2022

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

low_biom <- cbind(ceattle_biomass(run_wt001, "CEATTLE - 0.01% cannibalism"),
                  ceattle_biomass(run_wt005, "CEATTLE - 0.05% cannibalism"),
                  ceattle_biomass(run_wt01, "CEATTLE - 0.1% cannibalism"),
                  ceattle_biomass(run_wt05, "CEATTLE - 0.5% cannibalism"))
low_biom <- low_biom[, c(1:3, 6, 9, 12)]

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

plot_biom <- function(df) {
  wide <- cbind(df, nodiet_biom[, 3], ss_biom_all[, 3])
  colnames(wide)[(ncol(wide)-1):ncol(wide)] <- c("CEATTLE - no diet", "Stock Synthesis")
  biom <- melt(wide, id.vars = c("year", "type"))
  
  plot <- ggplot(biom, aes(x=year, y=value)) +
    geom_line(aes(color=variable)) +
    scale_color_viridis(discrete = TRUE, direction = -1) +  
    theme_sleek() +
    ylab("Biomass (mt)") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 1)
  
  return(plot)
}

# Low-cannibalism plot
low_biom_plot <- plot_biom(low_biom)
low_biom_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/low_intrasp_biomass.png", 
       low_biom_plot, width=200, height=150, units="mm", dpi=300)


# Plot recruitment ------------------------------------------------------------
nodiet_R <- read.csv("data/ceattle_nodiet_R.csv")
ss_R <- read.table("data/assessment/recruitment.txt")[15:57,]

low_R_data <- cbind(c(run_wt001$quantities$R), c(run_wt005$quantities$R),
                    c(run_wt01$quantities$R), c(run_wt05$quantities$R))
low_R_wide <- as.data.frame(cbind(years, low_R_data, nodiet_R))
colnames(low_R_wide) <- c("year",
                          "CEATTLE - 0.01% cannibalism",
                          "CEATTLE - 0.05% cannibalism",
                          "CEATTLE - 0.1% cannibalism",
                          "CEATTLE - 0.5% cannibalism",
                          "CEATTLE - no diet")
low_R <- melt(low_R_wide, id.vars = "year")

# Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS is 0)
ss_1 <- cbind(1981:2022, rep("SS + 1", (length(years)-1)), ss_R[1:42, 2])
colnames(ss_1) <- c("year", "variable", "value")

plot_R <- function(df) {
  df <- rbind(df, ss_1)
  df$value <- as.numeric(df$value)
  df$year <- as.numeric(df$year)
  
  plot <- ggplot(df, aes(x=year, y=value)) +
    geom_line(aes(color=variable)) +
    scale_color_viridis(discrete = TRUE, direction = -1) + 
    theme_sleek() +
    ylab("Recruitment") +
    labs(color = "model")
  
  return(plot)
}

low_R_plot <- plot_R(low_R)
low_R_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/low_intrasp_R.png", 
       low_R_plot, width=200, height=100, units="mm", dpi=300)


# Calculate and plot difference btw no diet & each cannibalism run ------------
nodiet_biom4 <- cbind(nodiet_biom[, 3], nodiet_biom[, 3],
                      nodiet_biom[, 3], nodiet_biom[, 3])
delta_biom_wide <- low_biom[44:86, c(3:6)] - nodiet_biom4
delta_biom_wide <- cbind(years, delta_biom_wide)
delta_biom <- melt(delta_biom_wide, id.vars = "years")

biom_difference <- ggplot(delta_biom2, aes(x=years, y=value)) +
  geom_line(aes(color=variable)) +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.336) +  # colors to match biomass plot  
  theme_sleek() +
  ylab("Biomass (mt)") + xlab("year")
  labs(color = "model") 
biom_difference

ggsave(filename="plots/CEATTLE/intraspecies predation/low_intrasp_biom_difference.png", 
       biom_difference, width=200, height=100, units="mm", dpi=300)


# Numbers-at-age for each model run -------------------------------------------
# Read in data from no diet CEATTLE run
nbyage_nodiet <- read.csv("data/ceattle_nodiet_nbyage.csv")
nbyage_nodiet <- cbind(nbyage_nodiet, rep("CEATTLE - no diet", nrow(nbyage_nodiet)))
colnames(nbyage_nodiet)[4] <- "model"

extract_nbyage <- function(run, name) {
  df <- as.data.frame(as.table(run$quantities$NByage))
  
  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(1980:2022)
  colnames(df) <- c("age", "year", "numbers")
  
  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"
  
  return(df)
}

nbyage_all <- rbind(extract_nbyage(run_wt001, "CEATTLE - 0.01% cannibalism"),
                    extract_nbyage(run_wt005, "CEATTLE - 0.05% cannibalism"),
                    extract_nbyage(run_wt01, "CEATTLE - 0.1% cannibalism"),
                    extract_nbyage(run_wt05, "CEATTLE - 0.5% cannibalism"),
                    nbyage_nodiet)

# Calculate mean numbers at age & plot
nbyage_mean <- nbyage_all %>% group_by(age, model) %>%
  summarize(mean_number = mean(numbers))

nbyage_plot_mean <- ggplot(nbyage_mean, aes(x=age, y=mean_number, fill=model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.166) +
  xlab("age") + ylab("numbers") 
nbyage_plot_mean

ggsave(filename = "plots/CEATTLE/intraspecies predation/nbyage_intrasp.png", 
       nbyage_plot_mean, width=200, height=120, units="mm", dpi=300)


# Run everything with higher amounts of cannibalism (not converging 100%) -----
# wt10 <- c(0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# wt30 <- c(0.0, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)
# wt50 <- c(0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
# wt80 <- c(0.0, 0.16, 0.24, 0.32, 0.40, 0.48, 0.56, 0.64, 0.72, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
# 
# run_wt05 <- run_ceattle(wt05, hake_intrasp)
# run_wt10 <- run_ceattle(wt10, hake_intrasp)
# run_wt30 <- run_ceattle(wt30, hake_intrasp)
# run_wt50 <- run_ceattle(wt50, hake_intrasp)
# run_wt80 <- run_ceattle(wt80, hake_intrasp)
# 
# high_test <- cbind(ceattle_biomass(run_wt05, "CEATTLE - 0.5% cannibalism"),
#                    ceattle_biomass(run_wt10, "CEATTLE - 10% cannibalism"),
#                    ceattle_biomass(run_wt30, "CEATTLE - 30% cannibalism"),
#                    ceattle_biomass(run_wt50, "CEATTLE - 50% cannibalism"),
#                    ceattle_biomass(run_wt80, "CEATTLE - 80% cannibalism"))
# high_test <- high_test[, c(1:3, 6, 9, 12, 15)]
# 
# 
# # High-cannibalism biomass plot
# high_biom_plot <- plot_biom(high_test)
# high_biom_plot
# 
# ggsave(filename="plots/CEATTLE/high_intrasp_biomass.png", high_biom_plot,
#        width=200, height=150, units="mm", dpi=300)
# 
# 
# # High-cannibalism recruitment plot
# high_R_data <- cbind(c(run_wt05$quantities$R), c(run_wt10$quantities$R),
#                      c(run_wt30$quantities$R), c(run_wt50$quantities$R),
#                      c(run_wt80$quantities$R))
# high_R_wide <- as.data.frame(cbind(years, high_R_data, nodiet_R))
# colnames(high_R_wide) <- c("year",
#                            "CEATTLE - 0.5% cannibalism",
#                            "CEATTLE - 10% cannibalism",
#                            "CEATTLE - 30% cannibalism",
#                            "CEATTLE - 50% cannibalism",
#                            "CEATTLE - 80% cannibalism",
#                            "CEATTLE - no diet")
# high_R <- melt(high_R_wide, id.vars = "year")
# 
# high_R_plot <- plot_R(high_R)
# high_R_plot
# 
# ggsave(filename="plots/CEATTLE/high_intrasp_R.png", high_R_plot,
#        width=200, height=100, units="mm", dpi=300)

