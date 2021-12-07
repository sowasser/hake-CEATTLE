# Script for constructing the bottom temperature timeseries for CEATTLE.
# Currently data are only available for survey years & are recorded as temp at
# 100m. Temperature was collected multiple times per year (probably each 
# survey station), but CEATTLE takes in a single data point per year.

library(dplyr)
library(ggplot2)

temp_all <- read.csv("data/temperature/temp_100_sophia.csv")

temp <- temp_all %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100))

temp_plot <- ggplot(temp, aes(x=year, y=mean_temp)) +
  geom_point() +
  theme_sleek()
# temp_plot

ggsave(filename="plots/survey_mean_temp.pdf", temp_plot,
       width=150, height=100, units="mm", dpi=300)