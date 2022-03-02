# Script for constructing the bottom temperature timeseries for CEATTLE.
# Currently data are only available for survey years & are recorded as temp at
# 100m. Temperature was collected multiple times per year (probably each 
# survey station), but CEATTLE takes in a single data point per year.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

temp_all <- read.csv("data/temperature/temp_100_sophia.csv")

# Temp summary stats ----------------------------------------------------------
overall_mean <- mean(temp_all$temp_100)  # 7.916°C
overall_median <- median(temp_all$temp_100)  # 7.763°C

overall_max <- max(temp_all$temp_100)  # 14.481°C
overall_min <- min(temp_all$temp_100)  # 5.734°C
overall_percentile <- quantile(temp_all$temp_100, 0.05, .95)  # 95% < 9.665°C


# Mean temperature per year for CEATTLE ---------------------------------------
# Including overall mean for missing years
temp_mean <- temp_all %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100))

missing_years <- cbind(c(1996, 1999, 1999, 2000, 2002, 2002, 2004, 2006, 2008,
                         2010, 2014, 2016, 2018, 2020), 
                       rep(overall_mean, length=14),
                       rep("overall mean", length=14))
colnames(missing_years) <- c("year", "temp", "data_type")

temp2 <- cbind(temp_mean, rep("year mean", length=13))
colnames(temp2) <- c("year", "temp", "data_type")

all_years <- rbind(temp2, missing_years)
all_years <- all_years[order(all_years$year), ]

all_years$temp <- as.numeric(all_years$temp)
all_years$year <- as.integer(all_years$year)

# Plot all mean temperatures 
mean_temp_plot <- ggplot(all_years, aes(x=year, y=temp)) +
  geom_line(color="gray") +
  geom_point(aes(color=data_type)) +
  scale_color_viridis(discrete = TRUE, direction=-1) +  # invert colors
  theme_sleek() +
  ylab("temperature")
# mean_temp_plot

ggsave(filename="plots/temperature/survey_mean_temp.png", mean_temp_plot,
       width=150, height=100, units="mm", dpi=300)

# Temp distribution & max for looking at hake thermal maximum -----------------
# Max temperature for each year
year_max <- temp_all %>% group_by(year) %>%
  filter(temp_100 == max(temp_100))

max_temp <- ggplot(year_max, aes(x=year, y=temp_100)) +
  geom_point() +
  theme_sleek() +
  ylab("temperature")
max_temp

ggsave(filename="plots/temperature/survey_max_temp.png", max_temp,
       width=150, height=100, units="mm", dpi=300)

# Temperature distribution overall and per year
temp_dist_all <- ggplot(temp_all, aes(x=temp_100)) +
  geom_histogram() +
  theme_sleek() +
  xlab("temperature") + ylab(" ")
temp_dist_all

ggsave(filename="plots/temperature/temp_dist_all.png", temp_dist_all,
       width=150, height=100, units="mm", dpi=300)

temp_dist_years <- ggplot(temp_all, aes(x=temp_100)) +
  geom_histogram() +
  theme_sleek() +
  xlab("temperature") + ylab(" ")
  facet_wrap(~year)
# temp_dist_years

ggsave(filename="plots/temperature/temp_dist_years.png", temp_dist_years,
       width=250, height=150, units="mm", dpi=300)


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

# Plot histogram of kriged temp values
temp_hake_hist <- ggplot(temp_hake, aes(x=temp_100_kriged)) +
  geom_histogram() +
  theme_sleek() +
  xlab("kriged temperature, biomass > 0") + ylab(" ")

ggsave(filename="plots/temperature/temp_hake_hist.png", temp_hake_hist,
       width=150, height=100, units="mm", dpi=300)


# Compare mean survey temp, mean kriged temp, kriged + hake biomass > 0 -------
temp_kriged_mean <- temp_kriged %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100_kriged))

temp_hake_mean <- temp_hake %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100_kriged))

temp_weighted <- temp_hake %>% group_by(year) %>%
  summarise(mean_temp = weighted.mean(temp_100_kriged, hake_biomass))

# Combine into 1 dataset with labeled data sources, then plot
means <- rbind(temp_mean, temp_kriged_mean, temp_hake_mean, temp_weighted)
means <- cbind(means, c(rep("survey", 13), 
                        rep("kriged, all", 12), 
                        rep("kriged, biomass > 0", 12),
                        rep("kriged, weighted", 12)))
colnames(means)[3] <- "dataset"

mean_temp_compared <- ggplot(means, aes(x=year, y=mean_temp)) +
  geom_point(aes(color=dataset), size=2) +
  geom_line(aes(color=dataset), size=1) +
  scale_color_viridis(discrete = TRUE) +   
  theme_sleek() +
  ylab("mean temperature")
mean_temp_compared

ggsave(filename="plots/temperature/mean_temp_compared.png", mean_temp_compared,
       width=160, height=100, units="mm", dpi=300)
