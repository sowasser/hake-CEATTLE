# Script for constructing the bottom temperature timeseries for CEATTLE.
# Currently data are only available for survey years & are recorded as temp at
# 100m. Temperature was collected multiple times per year (probably each 
# survey station), but CEATTLE takes in a single data point per year.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

survey_temp <- read.csv("data/temperature/temp_100_sophia.csv")[, -c(2:4)]
summer_ROMS <- read.csv("data/temperature/ROMS_summer_mean.csv")

missing_years <- c(1996, 1999, 1999, 2000, 2002, 2002, 2004, 2006, 2008,2010, 
                   2014, 2016, 2018, 2020)

# Combine summer ROMS mean with survey mean -----------------------------------
# Only keep years that don't overlap with the survey
# ROMS_nonsurvey <- filter(ROMS_mean, Year %in% c(1980:1994))

# Find mean from survey
survey_mean <- survey_temp %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100))

survey <- cbind(survey_mean, rep("survey", length(survey_mean$mean_temp)))
colnames(survey) <- c("year", "temp", "source")

ROMS <- cbind(summer_ROMS, rep("ROMS", length(summer_ROMS$mean_temp)))
colnames(ROMS) <- c("year", "mean_temp", "source")

# Combine together and sort by year
CEATTLE_temp <- rbind(ROMS, survey)
CEATTLE_temp <- CEATTLE_temp[order(CEATTLE_temp$year), ]

# Plot all mean temperatures 
mean_temp_plot <- ggplot(CEATTLE_temp, aes(x=year, y=mean_temp, color=source)) +
  geom_line(linetype="dotted") +
  geom_point() +
  scale_color_viridis(discrete = TRUE, direction=-1, begin=0.1, end=0.9) +  # invert colors
  theme_sleek() +
  ylab("temperature")
mean_temp_plot

ggsave(filename="plots/temperature/mean_temp.png", mean_temp_plot,
       width=150, height=100, units="mm", dpi=300)

# Temperature distribution per year -------------------------------------------
temp_dist_years <- ggplot(survey_temp, aes(x=temp_100)) +
  geom_histogram() +
  theme_sleek() +
  xlab("temperature") + ylab(" ") +
  facet_wrap(~year)
# temp_dist_years

ggsave(filename="plots/temperature/temp_dist_years.png", temp_dist_years,
       width=250, height=200, units="mm", dpi=300)


# Temperature where hake were found -------------------------------------------
# Temp here has been kriged to create a spatial surface for each year and then
# temp from this kriged grid was assigned to each hake biomass estimate.
# Info here: https://www.int-res.com/abstracts/meps/v639/p185-197/
temp_kriged <- read.csv("data/temperature/temp_100_matched_sophia.csv")

# Only include rows where hake were found
temp_hake <- temp_kriged %>%
  filter(hake_biomass > 0)

temp_hake_mean <- mean(temp_hake$temp_100_kriged)
temp_hake_max <- max(temp_hake$temp_100_kriged)
temp_hake_min <- min(temp_hake$temp_100_kriged)

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
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("temperature (°C)") + ylab(" ") 
temp_hake_hist

ggsave(filename="plots/temperature/temp_hake_hist.png", temp_hake_hist,
       width=160, height=100, units="mm", dpi=300)


# Compare mean survey temp, mean kriged temp, kriged + hake biomass > 0 -------
temp_kriged_mean <- temp_kriged %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100_kriged))

temp_hake_mean <- temp_hake %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100_kriged))

temp_weighted <- temp_hake %>% group_by(year) %>%
  summarise(mean_temp = weighted.mean(temp_100_kriged, hake_biomass))

# Combine into 1 dataset with labeled data sources, then plot
means <- rbind(ROMS[, 1:2], survey_mean, temp_kriged_mean, temp_weighted)
means <- cbind(means, c(rep("summer ROMS", 41),
                        rep("survey", 13), 
                        rep("kriged, all", 12), 
                        rep("kriged, biomass weighted", 12)))
colnames(means)[3] <- "dataset"

mean_temp_compared <- ggplot(means, aes(x=year, y=mean_temp)) +
  geom_point(aes(color=dataset), size=2) +
  geom_line(aes(color=dataset), size=1, linetype="dotted") +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.9) +   
  theme_sleek() +
  ylab("mean temperature")
mean_temp_compared

ggsave(filename="plots/temperature/mean_temp_compared.png", mean_temp_compared,
       width=160, height=100, units="mm", dpi=300)
