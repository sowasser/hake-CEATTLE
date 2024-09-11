# Script for constructing the bottom temperature timeseries for CEATTLE.
# Currently data are only available for survey years & are recorded as temp at
# 100m. Temperature was collected multiple times per year (probably each 
# survey station), but CEATTLE takes in a single data point per year.

library(dplyr)
library(ggplot2)
library(viridis)
library(ggsidekick)
library(here)
theme_set(theme_sleek())

TempDir <- here("data", "temperature")  # Directory for data

# Read in temperature datasets
survey_temp <- read.csv(here(TempDir, "temp_100_sophia.csv"))[, -c(2:4)]
summer_ROMS <- read.csv(here(TempDir, "ROMS_summer_mean.csv"))
ROMS <- read.csv(here(TempDir, "ROMS_mean.csv"))

survey_new <- readRDS(here(TempDir, "temperature_index.rds"))

min(ROMS$mean_temp)
max(ROMS$mean_temp)

missing_years <- c(1996, 1997, 1999, 1999, 2000, 2002, 2002, 2004, 2006, 2008,2010, 
                   2014, 2016, 2018, 2020)

# Combine summer ROMS mean with survey mean -----------------------------------
# Only keep years that don't overlap with the survey
# ROMS_nonsurvey <- filter(ROMS_mean, Year %in% c(1980:1994))

# Find mean from survey
survey_mean <- survey_temp %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100))

survey <- cbind(survey_mean, rep("survey", length(survey_mean$mean_temp)))
colnames(survey) <- c("year", "mean_temp", "source")

summer_ROMS <- cbind(summer_ROMS, rep("summer ROMS", length(summer_ROMS$mean_temp)))
colnames(summer_ROMS) <- c("year", "mean_temp", "source")

ROMS <- cbind(ROMS, rep("ROMS", length(ROMS$mean_temp)))
colnames(ROMS) <- c("year", "mean_temp", "source")

# Find difference between ROMS and survey
ROMS_surveyyears <- ROMS %>% 
  filter(!(year %in% missing_years) & year >= 1995) %>%
  mutate(ROMS_survey = survey_mean$mean_temp - mean_temp)
mean(ROMS_surveyyears$ROMS_survey)

# Get summer mean of new temp dataset 
survey_new_summer <- survey_new %>% 
  filter(month %in% 6:8) %>%
  group_by(year) %>%
  summarize(mean_temp = mean(mean_thetao)) %>%
  mutate(source = "Survey Area Mean")

write.csv(survey_new_summer, file = here(TempDir, "survey_summer_mean.csv"), row.names = FALSE)

# Plot just summer survey area mean
ggplot(survey_new_summer, aes(x = year, y = mean_temp)) +
  geom_point() +
  geom_line() +
  ylab("mean temperature")

# Combine together and sort by year
CEATTLE_temp <- rbind(ROMS, survey, survey_new_summer)
CEATTLE_temp <- CEATTLE_temp[order(CEATTLE_temp$year), ]

# Plot all mean temperatures 
mean_temp_plot <- ggplot(CEATTLE_temp, aes(x=year, y=mean_temp, color=source, shape = source)) +
  geom_point(size=2) +
  geom_line(linewidth=1, alpha = 0.3) +
  ylim(0, NA) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.6) +  # invert colors
  ylab("temperature")
mean_temp_plot

ggsave(filename="plots/temperature/mean_temp.png", mean_temp_plot,
       bg = "transparent", width=150, height=90, units="mm", dpi=300)

# Temperature distribution per year -------------------------------------------
temp_dist_years <- ggplot(survey_temp, aes(x=temp_100)) +
  geom_histogram() +
  xlab("temperature") + ylab(" ") +
  facet_wrap(~year)
# temp_dist_years

ggsave(filename="plots/temperature/temp_dist_years.png", temp_dist_years,
       bg = "transparent", width=250, height=200, units="mm", dpi=300)


# Temperature where hake were found -------------------------------------------
# Temp here has been kriged to create a spatial surface for each year and then
# temp from this kriged grid was assigned to each hake biomass estimate.
# Info here: https://www.int-res.com/abstracts/meps/v639/p185-197/
temp_kriged <- read.csv("data/temperature/temp_100_matched_sophia.csv")

# Only include rows where hake were found
temp_hake <- temp_kriged %>%
  filter(hake_biomass > 0)

# Mean temp for hake presence
mean(temp_hake$temp_100_kriged)
max(temp_hake$temp_100_kriged)
min(temp_hake$temp_100_kriged)

# Mean, weighted by hake biomass
weighted_mean <- weighted.mean(temp_hake$temp_100_kriged, temp_hake$hake_biomass)

# Plot histogram of kriged temp values vs. overall temperature 
temp_comp <- rbind(cbind(survey_temp$temp_100, 
                         rep("survey", length(survey_temp$temp_100))),
                   cbind(temp_hake$temp_100_kriged, 
                         rep("kriged, biomass > 0", length(temp_hake$temp_100_kriged))))
temp_comp <- as.data.frame(temp_comp)
colnames(temp_comp) <- c("temp", "source")
temp_comp$temp <- as.numeric(temp_comp$temp)
  
temp_hake_hist <- ggplot(temp_comp, aes(x=temp, fill=source)) +
  geom_histogram() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("temperature (°C)") + ylab(" ") 
temp_hake_hist

ggsave(filename="plots/temperature/temp_hake_hist.png", temp_hake_hist,
       bg = "transparent", width=160, height=100, units="mm", dpi=300)


# Compare mean survey temp, mean kriged temp, kriged + hake biomass > 0 -------
temp_kriged_mean <- temp_kriged %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100_kriged))

temp_hake_mean <- temp_hake %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100_kriged))

temp_weighted <- temp_hake %>% group_by(year) %>%
  summarise(mean_temp = weighted.mean(temp_100_kriged, hake_biomass))

# Combine into 1 dataset with labeled data sources, then plot
means <- rbind(ROMS[, 1:2], survey_mean, temp_weighted)  # subset to model years
means <- cbind(means, c(rep("ROMS", length(1:41)),
                        rep("Survey", 13), 
                        rep("Kriged, biomass-weighted", 12)))
colnames(means)[3] <- "dataset"
means$dataset <- factor(means$dataset, levels = c("ROMS", "Survey", "Kriged, biomass-weighted"))

mean_temp_compared <- ggplot(means, aes(x=year, y=mean_temp)) +
  geom_point(aes(color=dataset), size=2) +
  geom_line(aes(color=dataset), linewidth=1, alpha = 0.3) +
  ylim(0, NA) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.9) +   
  ylab("Mean Temperature") + ylab("Year") + labs(color = "Dataset")
mean_temp_compared

ggsave(filename="plots/temperature/mean_temp_compared.png", mean_temp_compared,
       bg = "transparent", width=180, height=90, units="mm", dpi=300)
