# Run CEATTLE with intraspecies-predation proportions calculated from diet 
# database going back to 1980.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_220713.xlsx")

intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                                 inits = NULL, # Initial parameters = 0
                                 file = NULL, # Don't save
                                 # debug = 1, # 1 = estimate, 0 = don't estimate
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # Single species mode
                                 phase = "default")

# Check what all comes out of CEATTLE
ceattle_stuff <- intrasp_run$quantities

# # Rceattle diagmostics
# plot_index(intrasp_run)
# plot_catch(intrasp_run)
# plot_selectivity(intrasp_run)
# plot_mortality(intrasp_run)
# plot_indexresidual(intrasp_run)
# plot_logindex(intrasp_run)
# plot_recruitment(intrasp_run, add_ci = TRUE)
# plot_comp(intrasp_run)
# plot_srv_comp(intrasp_run)


### Plot biomass in comparison to no diet & assessment ------------------------
years <- 1980:2022
test <- intrasp_run$quantities$biomass
# Pull out SSB & overall biomass from CEATTLE runs
ceattle_biomass <- function(run, name) {
  ssb <- (c(run$quantities$biomassSSB) * 2)
  biom <- c(run$quantities$biomass)
  wide <- as.data.frame(cbind(years, ssb, biom))
  colnames(wide) <- c("year", "SSB", "Total Biomass")
  all_biom <- melt(wide, id.vars = "year")
  all_biom2 <- cbind(all_biom, 
                     error = rep(0, length(all_biom$year)),  # add error, 0 for now
                     model = rep(name, length(all_biom$year)))  
  colnames(all_biom2)[2:3] <- c("type", "value")
  
  return(all_biom2)
}

biomass <- ceattle_biomass(intrasp_run, "CEATTLE - cannibalism")
write.csv(biomass, "data/ceattle_intrasp_biomass.csv", row.names = FALSE)

# Read in other model runs for comparison and plot
plot_biomass <- function() {
  nodiet_biom <- read.csv("data/ceattle_nodiet_biom.csv")
  colnames(nodiet_biom)[3] <- "value"
  nodiet_biom <- cbind(nodiet_biom, 
                       error = rep(0, length(nodiet_biom$year)),  # add error, 0 for now
                       model = rep("CEATTLE - no diet", length(nodiet_biom$year)))  
  
  # Pull out SSB & total biomass from stock synthesis & combine, remove pre-1980
  ss3_ssb <- cbind(read.table("data/assessment/ssb.txt")[15:57, 2:3], type = rep("SSB", length(15:57)))
  ss3_biomass <- cbind(read.table("data/assessment/biomass.txt")[15:57, 2:3], type = rep("Total Biomass", length(15:57)))
  
  ss3_biom <- as.data.frame(cbind(year = rep(years, 2), 
                                  rbind(ss3_ssb, ss3_biomass),
                                  model = rep("Stock Synthesis", length(ss3_ssb$V2) * 2)))
  colnames(ss3_biom)[2:3] <- c("value", "error")
  ss3_biom <- ss3_biom[, c(1, 4, 2, 3, 5)]
  
  biom_all <- rbind(biomass, nodiet_biom, ss3_biom)
  
  biom_plot <- ggplot(biom_all, aes(x=year, y=value, color = model, fill = model)) +
    geom_line() +
    geom_ribbon(aes(ymin=(value-error), ymax=(value+error)), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    ylab("Biomass (mt)") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 1)
  
  return(biom_plot)
}

biom_plot <- plot_biomass()

biom_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/intrasp_biomass.png", biom_plot, 
       bg = "transparent", width=200, height=150, units="mm", dpi=300)


### Plot recruitment ----------------------------------------------------------
recruitment <- c(intrasp_run$quantities$R)
write.csv(recruitment, "data/ceattle_intrasp_R.csv", row.names = FALSE)

plot_R <- function() {
  nodiet_R <- read.csv("data/ceattle_nodiet_R.csv")
  ss3_R <- read.table("data/assessment/recruitment.txt")[15:57,]
  
  R_wide <- data.frame(year = years, recruitment, nodiet_R)
  colnames(R_wide)[2:3] <- c("CEATTLE - cannibalism", "CEATTLE - nodiet")
  
  R <- melt(R_wide, id.vars = "year")
  
  # Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS3 is 0)
  ss3_1 <- as.data.frame(cbind(year = 1981:2022, 
                               variable = rep("SS3 + 1", (length(1981:2022))), 
                               value = ss3_R[1:42, 2],
                               error = ss3_R[1:42, 3]))
  
  R_all <- rbind(cbind(R, error = rep(0, length(2 * R$value))), 
                 ss3_1)
  R_all$value <- as.numeric(R_all$value)
  R_all$year <- as.numeric(R_all$year)
  R_all$error <- as.numeric(R_all$error)
  
  R_plot <- ggplot(R_all, aes(x=year, y=value, color=variable, fill=variable)) +
    geom_line(aes(color=variable)) +
    geom_ribbon(aes(ymin=(value-error), ymax=(value+error)), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    ylab("Recruitment") +
    labs(color = "model") +labs(fill = "model")
  return(R_plot)
}

R_plot <- plot_R()
R_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/intrasp_R.png", R_plot, 
       bg = "transparent", width=200, height=100, units="mm", dpi=300)


### Numbers-at-age for each model run -----------------------------------------
# Extract numbers at age for intraspecies predation model run
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

nbyage <- extract_nbyage(intrasp_run, "CEATTLE - cannibalism")
write.csv(nbyage, "data/ceattle_intrasp_nbyage.csv", row.names = FALSE)

plot_nbyage <- function(output) {
  # Read in data from no diet CEATTLE run
  nbyage_nodiet <- read.csv("data/ceattle_nodiet_nbyage.csv")
  nbyage_nodiet <- cbind(nbyage_nodiet, rep("CEATTLE - no diet", nrow(nbyage_nodiet)))
  colnames(nbyage_nodiet)[4] <- "model"
  
  # Read in data from SS3 & average beginning & middle of the year
  nbyage_ss3_all <- read.csv("data/assessment/nbyage.csv")
  colnames(nbyage_ss3_all) <- c("year", "timing", c(0:20))
  
  nbyage_ss3_wide <- nbyage_ss3_all %>%
    group_by(year) %>%
    summarize_at(vars("0":"20"), mean)
  
  nbyage_ss3 <- melt(nbyage_ss3_wide[, -2], id.vars = "year")
  nbyage_ss3 <- cbind(nbyage_ss3, rep("Stock Synthesis", length(nbyage_ss3$year)))
  colnames(nbyage_ss3)[2:4] <- c("age", "numbers", "model")
  
  # Combine with nbyage from intrasp run
  nbyage_all <- rbind(nbyage, nbyage_nodiet, nbyage_ss3)
  
  # Set 15 as accumulation age
  nbyage_all$age[as.numeric(nbyage_all$age) > 15] <- 15
  
  # Calculate mean numbers at age & plot
  nbyage_mean <- nbyage_all %>% 
    group_by(age, model) %>%
    summarize(mean = mean(numbers), sd = sd(numbers)) 
  
  # Plot mean nbyage across years
  nbyage_plot_mean <- ggplot(nbyage_mean, aes(x=age, y=mean, fill=model, color=model)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) +  # only upper error bars
    scale_x_discrete(labels = c(1:14, "15+")) +
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    xlab("age") + ylab("numbers") 
  
  # Plot yearly nbyage
  nbyage_plot_yearly <- ggplot(nbyage_all, aes(x=age, y=numbers, fill=model)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete(labels = c(1:14, "15+")) +
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    xlab("age") + ylab("numbers") +
    facet_wrap(~ year)
  
  # Conditionally return plots
  if(output == "mean") {return(nbyage_plot_mean)}
  if(output == "yearly") {return(nbyage_plot_yearly)}
}

nbyage_plot_mean <- plot_nbyage(output = "mean")
nbyage_plot_mean

# nbyage_plot_yearly <- plot_nbyage(output = "yearly")
# nbyage_plot_yearly

ggsave(filename = "plots/CEATTLE/intraspecies predation/nbyage_intrasp.png", nbyage_plot_mean, 
       bg = "transparent", width=200, height=120, units="mm", dpi=300)


### Compare survey biomass estimate from CEATTLE to true values ---------------
intrasp_srv <- data.frame(year = 1995:2019,
                          biomass = intrasp_run$quantities$srv_bio_hat,
                          log_sd = intrasp_run$quantities$srv_log_sd_hat,
                          model = rep("CEATTLE - cannibalism", length(1995:2019)))

plot_survey <- function() {
  nodiet_srv <- read.csv("data/ceattle_nodiet_survey.csv")
  nodiet_srv <- cbind(nodiet_srv, model = rep("CEATTLE - no diet", length(nodiet_srv$year)))
  
  survey <- read.csv("data/assessment/survey_data.csv")
  survey <- cbind(survey, model = rep("assessment", length(survey$year)))
  
  survey_all <- rbind(intrasp_srv, nodiet_srv, survey)
  
  survey_plot <- ggplot(survey_all, aes(x=year, y=biomass, color=model)) +
    geom_line(alpha = 0.3) +
    geom_point() +
    # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    xlab("year") + ylab("survey biomass") 
  
  return(survey_plot)
}

survey_plot <- plot_survey()
survey_plot

ggsave(filename = "plots/CEATTLE/intraspecies predation/survey_biomass.png", survey_plot, 
       bg = "transparent", width=200, height=120, units="mm", dpi=300)


### Compare predation mortality (M2) ------------------------------------------
# M2 <- intrasp_run$quantities$M2
# M2_prop <- as.data.frame(intrasp_run$quantities$M2_prop)
# M <- as.data.frame(intrasp_run$quantities$M)
# 
# # We can plot all runs
# mod_list <- list(intrasp_run)
# mod_names <- c("Intraspecies")
# 
# # Plot biomass trajectory
# plot_biomass(Rceattle = mod_list, model_names = mod_names)
# plot_depletionSSB(Rceattle = mod_list, model_names = mod_names)
# plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
# 
# # Plot mortality and predation
# plot_b_eaten(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey
# plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names) # Biomass eaten as prey by each predator
plot_mortality(Rceattle = intrasp_run, type = 1) # Mortality-at-age time series
# 
# # Run diagnostics
# plot_selectivity(Rceattle = intrasp_run)
# plot_comp(intrasp_run) # Fitted survey composition data
# plot_index(intrasp_run) # Fitted indices of abundance
# plot_catch(intrasp_run) # Fitted catch series
