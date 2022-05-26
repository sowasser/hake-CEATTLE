# Run CEATTLE with intraspecies-predation proportions calculated from diet 
# database going back to 1980.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

hake_intrasp <- Rceattle::read_data( file = "data/hake_intrasp_220524.xlsx")

intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                                 inits = NULL, # Initial parameters = 0
                                 file = NULL, # Don't save
                                 # debug = 1, # 1 = estimate, 0 = don't estimate
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 0, # Single species mode
                                 phase = "default")

# Check what all comes out of CEATTLE
ceattle_stuff <- intrasp_run$quantities


### Plot biomass in comparison to no diet & assessment ------------------------
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

biomass <- ceattle_biomass(intrasp_run, "CEATTLE - cannibalism")

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
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    theme_sleek() +
    ylab("Biomass (mt)") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 1)
  
  return(plot)
}

# Low-cannibalism plot
biom_plot <- plot_biom(biomass)
biom_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/intrasp_biomass.png", 
       biom_plot, width=200, height=150, units="mm", dpi=300)


### Plot recruitment ----------------------------------------------------------
nodiet_R <- read.csv("data/ceattle_nodiet_R.csv")
ss3_R <- read.table("data/assessment/recruitment.txt")[15:57,]

R_wide <- data.frame(year = years, R_intrasp = c(intrasp_run$quantities$R), nodiet_R)
colnames(R_wide)[3] <- "R_nodiet"

R <- melt(R_wide, id.vars = "year")

# Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS is 0)
ss3_1 <- cbind(year = 1980:2022, 
              variable = rep("SS + 1", (length(years)-1)), 
              value = ss3_R[1:42, 2],
              error = ss3_R[1:42, 3])

plot_R <- function(df) {
  df <- rbind(cbind(df, error = rep(0, length(2 * df$value))), ss3_1)
  df$value <- as.numeric(df$value)
  df$year <- as.numeric(df$year)
  df$error <- as.numeric(df$error)
  
  plot <- ggplot(df, aes(x=year, y=value, color=variable, fill=variable)) +
    geom_line(aes(color=variable)) +
    geom_ribbon(aes(ymin=(value-error), ymax=(value+error)), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    theme_sleek() +
    ylab("Recruitment") +
    labs(color = "model") +labs(fill = "model")
  
  return(plot)
}

R_plot <- plot_R(R)
R_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/intrasp_R.png", 
       R_plot, width=200, height=100, units="mm", dpi=300)


### Numbers-at-age for each model run -----------------------------------------
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

nbyage_ss <- melt(nbyage_ss3_wide[, -2], id.vars = "year")
nbyage_ss <- cbind(nbyage_ss, rep("Stock Synthesis", length(nbyage_ss$year)))
colnames(nbyage_ss)[2:4] <- c("age", "numbers", "model")


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

nbyage_all <- rbind(extract_nbyage(intrasp_run, "CEATTLE - cannibalism"),
                    nbyage_nodiet, nbyage_ss)

# Set 15 as accumulation age
nbyage_all$age[as.numeric(nbyage_all$age) > 15] <- 15

# Calculate mean numbers at age & plot
nbyage_mean <- nbyage_all %>% 
  group_by(age, model) %>%
  summarize(mean = mean(numbers), sd = sd(numbers)) 

nbyage_plot_mean <- ggplot(nbyage_mean, aes(x=age, y=mean, fill=model, color=model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) +  # only upper error bars
  theme_sleek() +
  scale_x_discrete(labels = c(1:14, "15+")) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  xlab("age") + ylab("numbers") 
nbyage_plot_mean

ggsave(filename = "plots/CEATTLE/intraspecies predation/nbyage_intrasp.png", 
       nbyage_plot_mean, width=200, height=120, units="mm", dpi=300)


nbyage_plot_yearly <- ggplot(nbyage_all, aes(x=age, y=numbers, fill=model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_sleek() +
  scale_x_discrete(labels = c(1:14, "15+")) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  xlab("age") + ylab("numbers") +
  facet_wrap(~ year)
# nbyage_plot_yearly


### Compare survey biomass estimate from CEATTLE to true values ---------------
nodiet_srv <- read.csv("data/ceattle_nodiet_survey.csv")
nodiet_srv <- cbind(nodiet_srv, model = rep("CEATTLE - no diet", length(nodiet_srv$year)))

survey <- read.csv("data/assessment/survey_data.csv")
survey <- cbind(survey, model = rep("assessment", length(survey$year)))

intrasp_srv <- data.frame(year = 1995:2019,
                          biomass = intrasp_run$quantities$srv_bio_hat,
                          log_sd = intrasp_run$quantities$srv_log_sd_hat,
                          model = rep("CEATTLE - cannibalism", length(1995:2019)))

survey_all <- rbind(intrasp_srv, nodiet_srv, survey)

survey_plot <- ggplot(survey_all, aes(x=year, y=biomass, color=model)) +
  geom_line(linetype = "dotted") +
  geom_point() +
  # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
  theme_sleek() +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  xlab("year") + ylab("survey biomass") 
survey_plot

ggsave(filename = "plots/CEATTLE/intraspecies predation/survey_biomass.png", 
       survey_plot, width=200, height=120, units="mm", dpi=300)


