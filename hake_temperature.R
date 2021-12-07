# Script for constructing the bottom temperature timeseries for CEATTLE.
# Currently data are only available for survey years & are recorded as temp at
# 100m. Temperature was collected multiple times per year (probably each 
# survey station), but CEATTLE takes in a single data point per year.

library(dplyr)
library(ggplot2)

temp_all <- read.csv("data/temperature/temp_100_sophia.csv")

temp <- temp_all %>% group_by(year) %>%
  summarise(mean_temp = mean(temp_100))

# Look at missing years - from CEATTLE: Will use the mean for missing years.
all_mean <- mean(temp$mean_temp)

missing_years <- cbind(c(1996, 1999, 1999, 2000, 2002, 2002, 2004, 2006, 2008,
                         2010, 2014, 2016, 2018, 2020), 
                       rep(all_mean, length=14),
                       rep("overall mean", length=14))
colnames(missing_years) <- c("year", "temp", "data_type")

temp2 <- cbind(temp, rep("year mean", length=13))
colnames(temp2) <- c("year", "temp", "data_type")

all_years <- rbind(temp2, missing_years)
all_years <- all_years[order(all_years$year), ]

temp_plot <- ggplot(all_years, aes(x=year, y=temp, color=data_type)) +
  geom_point() +
  theme_sleek()
temp_plot

ggsave(filename="plots/survey_mean_temp.pdf", temp_plot,
       width=150, height=100, units="mm", dpi=300)