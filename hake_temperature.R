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
  theme_sleek()
# mean_temp_plot

ggsave(filename="plots/temperature/survey_mean_temp.pdf", mean_temp_plot,
       width=150, height=100, units="mm", dpi=300)

# Temp distribution & max for looking at hake thermal maximum -----------------
# Max temperature for each year
year_max <- temp_all %>% group_by(year) %>%
  filter(temp_100 == max(temp_100))

max_temp <- ggplot(year_max, aes(x=year, y=temp_100)) +
  geom_point() +
  theme_sleek()
max_temp

ggsave(filename="plots/temperature/survey_max_temp.pdf", max_temp,
       width=150, height=100, units="mm", dpi=300)

# Temperature distribution overall and per year
temp_dist_all <- ggplot(temp_all, aes(x=temp_100)) +
  geom_histogram() +
  theme_sleek() 
temp_dist_all

ggsave(filename="plots/temperature/temp_dist_all.pdf", temp_dist_all,
       width=150, height=100, units="mm", dpi=300)

temp_dist_years <- ggplot(temp_all, aes(x=temp_100)) +
  geom_histogram() +
  theme_sleek() +
  facet_wrap(~year)
# temp_dist_years

ggsave(filename="plots/temperature/temp_dist_years.pdf", temp_dist_years,
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

# Plot histogram of kriged temp values
temp_kriged <- ggplot(temp_hake, aes(x=temp_100_kriged)) +
  geom_histogram() +
  theme_sleek()

ggsave(filename="plots/temperature/temp_kriged.pdf", temp_kriged,
       width=150, height=100, units="mm", dpi=300)
