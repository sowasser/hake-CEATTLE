# Run CEATTLE for Pacific hake while testing different proportions of
# intraspecies predation

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

hake_intrasp <- Rceattle::read_data( file = "data/hake_intrasp_220524.xlsx")

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
wt01 <- c(0.0, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001)
wt05 <- c(0.0, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005)
wt1 <- c(0.0, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
wt5 <- c(0.0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
wt10 <- c(0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
wt30 <- c(0.0, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)
wt50 <- c(0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
wt80 <- c(0.0, 0.16, 0.24, 0.32, 0.40, 0.48, 0.56, 0.64, 0.72, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)

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
run_wt01 <- run_ceattle(wt01, hake_intrasp)
run_wt05 <- run_ceattle(wt05, hake_intrasp)
run_wt1 <- run_ceattle(wt1, hake_intrasp)
run_wt5 <- run_ceattle(wt5, hake_intrasp)
run_wt10 <- run_ceattle(wt10, hake_intrasp)
run_wt30 <- run_ceattle(wt30, hake_intrasp)
run_wt50 <- run_ceattle(wt50, hake_intrasp)
run_wt80 <- run_ceattle(wt80, hake_intrasp)

# Check what all comes out of CEATTLE
ceattle_stuff <- run_wt05$quantities


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

test_biom <- cbind(ceattle_biomass(run_wt05, "CEATTLE - 0.5% cannibalism"),
                   ceattle_biomass(run_wt10, "CEATTLE - 10% cannibalism"),
                   ceattle_biomass(run_wt50, "CEATTLE - 50% cannibalism"),
                   ceattle_biomass(run_wt80, "CEATTLE - 80% cannibalism"))
test_biom <- test_biom[, c(1:3, 6, 9, 12)]

# Read in no diet data
nodiet_biom <- read.csv("data/ceattle_nodiet_biom.csv")
colnames(nodiet_biom)[3] <- "CEATTLE - no diet"

# Pull out SSB & total biomass from stock synthesis & combine, remove pre-1980
ss3_ssb_werror <- read.table("data/assessment/ssb.txt")
ss3_ssb <- ss3_ssb_werror[15:57, 2]

ss3_biomass <- read.table("data/assessment/biomass.txt")
ss3_biom <- ss3_biomass[15:57, 2]

ss3_biom_wide <- as.data.frame(cbind(years, ss3_ssb, ss3_biom))
colnames(ss3_biom_wide) <- c("year", "SSB", "total biomass")
ss3_biom_all <- melt(ss3_biom_wide, id.vars = "year")

plot_biom <- function(df) {
  wide <- cbind(df, nodiet_biom[, 3], ss3_biom_all[, 3])
  colnames(wide)[(ncol(wide)-1):ncol(wide)] <- c("CEATTLE - no diet", "Stock Synthesis")
  biom <- melt(wide, id.vars = c("year", "type"))
  
  plot <- ggplot(biom, aes(x=year, y=value)) +
    geom_line(aes(color=variable)) +
    scale_color_viridis(discrete = TRUE, direction = -1) +  
    ylab("Biomass (mt)") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 1)
  
  return(plot)
}

# Low-cannibalism plot
test_biom_plot <- plot_biom(test_biom)
test_biom_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/test_intrasp_biomass.png", test_biom_plot, 
       bg = "transparent", width=200, height=150, units="mm", dpi=300)


# Plot recruitment ------------------------------------------------------------
nodiet_R <- read.csv("data/ceattle_nodiet_R.csv")
ss3_R <- read.table("data/assessment/recruitment.txt")[15:57,]

R_test_all <- cbind(c(run_wt05$quantities$R), c(run_wt10$quantities$R),  
                    c(run_wt50$quantities$R), c(run_wt80$quantities$R))
R_test_wide <- as.data.frame(cbind(years, R_test_all, nodiet_R))
colnames(R_test_wide) <- c("year",
                           "CEATTLE - 0.5% cannibalism",
                           "CEATTLE - 10% cannibalism",
                           "CEATTLE - 50% cannibalism",
                           "CEATTLE - 80% cannibalism",
                           "CEATTLE - no diet")
R_test <- melt(R_test_wide, id.vars = "year")

# Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS is 0)
ss3_1 <- cbind(1981:2022, rep("SS3 + 1", (length(years)-1)), ss3_R[1:42, 2])
colnames(ss3_1) <- c("year", "variable", "value")

plot_R <- function(df) {
  df <- rbind(df, ss3_1)
  df$value <- as.numeric(df$value)
  df$year <- as.numeric(df$year)
  
  plot <- ggplot(df, aes(x=year, y=value)) +
    geom_line(aes(color=variable)) +
    scale_color_viridis(discrete = TRUE, direction = -1) + 
    ylab("Recruitment") +
    labs(color = "model")
  
  return(plot)
}

test_R_plot <- plot_R(R_test)
test_R_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/test_intrasp_R.png", test_R_plot, 
       bg = "transparent", width=200, height=100, units="mm", dpi=300)


# # Calculate and plot difference btw no diet & each cannibalism run ------------
# nodiet_biom4 <- cbind(nodiet_biom[, 3], nodiet_biom[, 3],
#                       nodiet_biom[, 3], nodiet_biom[, 3])
# delta_biom_wide <- low_biom[44:86, c(3:6)] - nodiet_biom4
# delta_biom_wide <- cbind(years, delta_biom_wide)
# delta_biom <- melt(delta_biom_wide, id.vars = "years")
# 
# biom_difference <- ggplot(delta_biom2, aes(x=years, y=value)) +
#   geom_line(aes(color=variable)) +
#   scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.336) +  # colors to match biomass plot  
#   theme_sleek() +
#   ylab("Biomass (mt)") + xlab("year")
#   labs(color = "model") 
# biom_difference
# 
# ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/low_intrasp_biom_difference.png", 
#        biom_difference, width=200, height=100, units="mm", dpi=300)


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

nbyage_test_all <- rbind(extract_nbyage(run_wt05, "CEATTLE - 0.5% cannibalism"),
                         extract_nbyage(run_wt10, "CEATTLE - 10% cannibalism"),
                         extract_nbyage(run_wt50, "CEATTLE - 50% cannibalism"),
                         extract_nbyage(run_wt80, "CEATTLE - 80% cannibalism"),
                         nbyage_nodiet)

# Set 15 as accumulation age
nbyage_test_all$age[as.numeric(nbyage_test_all$age) > 15] <- 15

# Calculate mean numbers at age & plot
nbyage_test_mean <- nbyage_test_all %>% group_by(age, model) %>%
  summarize(mean_number = mean(numbers))

test_nbyage_plot <- ggplot(nbyage_test_mean, aes(x=age, y=mean_number, fill=model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(labels = c(1:14, "15+")) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.166) +
  xlab("age") + ylab("numbers") 
test_nbyage_plot

ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_intrasp_nbyage.png", test_nbyage_plot, 
       bg = "transparent", width=200, height=120, units="mm", dpi=300)


### Compare survey biomass estimate from CEATTLE to true values ---------------
nodiet_srv <- read.csv("data/ceattle_nodiet_survey.csv")
nodiet_srv <- cbind(nodiet_srv, model = rep("CEATTLE - no diet", length(nodiet_srv$year)))

survey <- read.csv("data/assessment/survey_data.csv")
survey <- cbind(survey, model = rep("SS3", length(survey$year)))

extract_srv <- function(run, name){
  df <- data.frame(year = 1995:2019,
                   biomass = run$quantities$srv_bio_hat,
                   log_sd = run$quantities$srv_log_sd_hat,
                   model = rep(name, length(1995:2019)))
  return(df)
}

srv_test <- rbind(extract_srv(run_wt05, "CEATTLE - 0.5% cannibalism"),
                  extract_srv(run_wt10, "CEATTLE - 10% cannibalism"),
                  extract_srv(run_wt50, "CEATTLE - 50% cannibalism"),
                  extract_srv(run_wt80, "CEATTLE - 80% cannibalism"),
                  nodiet_srv,
                  survey)

test_survey_plot <- ggplot(srv_test, aes(x=year, y=biomass, color=model)) +
  geom_line(linetype = "dotted") +
  geom_point() +
  # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  xlab("year") + ylab("survey biomass") 
test_survey_plot

ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_intrasp_survey.png", test_survey_plot, 
       bg = "transparent", width=200, height=120, units="mm", dpi=300)
