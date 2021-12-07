# Script for constructing the bottom temperature timeseries for CEATTLE.
# Currently data are only available for survey years & are recorded as temp at
# 100m. Temperature was collected multiple times per year (probably each 
# survey station), but CEATTLE takes in a single data point per year.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

temp_all <- read.csv("data/temperature/temp_100_sophia.csv")

# Mean temperature per year for CEATTLE ---------------------------------------
# Iincluding overall mean for missing years
temp_mean <- temp_all %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100))

all_mean <- mean(temp_mean$mean_temp)

missing_years <- cbind(c(1996, 1999, 1999, 2000, 2002, 2002, 2004, 2006, 2008,
                         2010, 2014, 2016, 2018, 2020), 
                       rep(all_mean, length=14),
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

ggsave(filename="plots/survey_mean_temp.pdf", mean_temp_plot,
       width=150, height=100, units="mm", dpi=300)

# Temp distribution & max for looking at hake thermal maximum -----------------
# overall maximum
overall_max <- max(temp_all$temp_100)  # 14.481°C

# Max temperature for each year
year_max <- temp_all %>% group_by(year) %>%
  filter(temp_100 == max(temp_100))

max_temp <- ggplot(year_max, aes(x=year, y=temp_100)) +
  geom_point() +
  theme_sleek()
max_temp

ggsave(filename="plots/survey_max_temp.pdf", max_temp,
       width=150, height=100, units="mm", dpi=300)

# Temperature distribution overall and per year
temp_dist_all <- ggplot(temp_all, aes(x=temp_100)) +
  geom_histogram() +
  theme_sleek() 
temp_dist_all

ggsave(filename="plots/temp_dist_all.pdf", temp_dist_all,
       width=150, height=100, units="mm", dpi=300)

temp_dist_years <- ggplot(temp_all, aes(x=temp_100)) +
  geom_histogram() +
  theme_sleek() +
  facet_wrap(~year)
# temp_dist_years

ggsave(filename="plots/temp_dist_years.pdf", temp_dist_years,
       width=150, height=100, units="mm", dpi=300)
