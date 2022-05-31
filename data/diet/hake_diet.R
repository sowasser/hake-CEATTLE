# Script for analyzing diet data for hake. First is for visualizing diet 
# proportion by observation and by weight for input in CEATTLE. Second is 
# analyzing the hake diet data from the US acoustic trawl survey. Data are 
# split between the predators (all hake) and the prey that are hake.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(fishmethods)
library(tidyr)
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(rgeos)


# Acoustic trawl survey data --------------------------------------------------
stomach_summary <- read.csv("data/diet/old/hake_stomachs_summary.csv")
full_stomachs <- read.csv("data/diet/old/full_stomachs.csv")
hake_predator <- read.csv("data/diet/old/hake_predator_all.csv")
all_stomachs <- read.csv("data/diet/old/hake_stomachs_all.csv")

# Total stomachs per year -----------------------------------------------------
per_year <- ggplot(stomach_summary, aes(x=year, y=num_stomachs_tow, fill=type, color=type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("year") + ylab("total stomachs")
per_year

ggsave(filename="plots/diet/old/stomachs_per_year.png", per_year,
       width=150, height=100, units="mm", dpi=300)

# Full & empty stomach locations ----------------------------------------------
hake_locations <- hake_predator %>%
  group_by(year, td_latitude, td_longitude, common_name) %>%
  summarize(n = n()) %>%
  arrange(common_name)

colnames(hake_locations)[4] <- "contents"

# Create a plot of location of observations by latitude and longitude
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry

hake_locations_yearly <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = hake_locations, aes(x = td_longitude, y = td_latitude, size = n, color = contents)) +
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) +
  scale_x_continuous(breaks = seq(-135, -120, by = 10)) +
  scale_y_continuous(breaks = seq(35, 55, by = 10)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab(" ") + ylab(" ") +
  facet_wrap(~year, ncol = 5)
hake_locations_yearly

ggsave(filename = "plots/diet/old/cannibalism_locations.png", hake_locations_yearly, 
       width=200, height=150, units="mm", dpi=300)


# Full stomach pred ages ------------------------------------------------------
full_pred_age <- ggplot(full_stomachs, aes(x=age)) +
  geom_histogram() +
  theme_sleek() +
  xlab("predator age") + ylab(" ")
full_pred_age

ggsave(filename="plots/diet/old/full_predator_age.png", full_pred_age,
       width=150, height=100, units="mm", dpi=300)


# Empty vs. hake_containing stomachs ------------------------------------------
empty_hake_length <- ggplot(hake_predator, aes(x=fork_length.cm, fill=scientific_name)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, labels = c("empty", "hake"), begin = 0.1, end = 0.9) +
  xlab("fork length (cm)") + ylab(" ") +
  labs(fill = "contents")
empty_hake_length

ggsave(filename="plots/diet/old/emptyvshake_length.png", empty_hake_length,
       width=150, height=100, units="mm", dpi=300)

empty_hake_age <- ggplot(hake_predator, aes(x=age, fill=scientific_name)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, labels = c("empty", "hake"), begin = 0.1, end = 0.9) +
  xlab("predator age") + ylab(" ") +
  labs(fill = "contents")
empty_hake_age

ggsave(filename="plots/diet/old/emptyvshake_age.png", empty_hake_age,
       width=150, height=100, units="mm", dpi=300)
  

# Predator hake dynamics ------------------------------------------------------
pred_lengths <- ggplot(hake_predator, aes(x=fork_length.cm, fill=sex)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("fork length (cm)") + ylab(" ") 
pred_lengths

ggsave(filename="plots/diet/old/all_pred_lengths.png", pred_lengths,
       width=150, height=100, units="mm", dpi=300)

pred_age <- ggplot(hake_predator, aes(x=age, fill=sex)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("age") + ylab(" ")
pred_age

ggsave(filename="plots/diet/old/pred_age.png", pred_age,
       width=150, height=100, units="mm", dpi=300)

pred_sex <- ggplot(hake_predator, aes(x=year, fill=sex)) +
  geom_bar(stat="count") +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("year") + ylab(" ")
pred_sex

ggsave(filename="plots/diet/old/pred_sex.png", pred_sex,
       width=150, height=100, units="mm", dpi=300)


# Length & weight of predator/prey --------------------------------------------
lengths <- as.data.frame(
  rbind(cbind(na.omit(full_stomachs$content_length.cm), 
              rep("prey", times = length(na.omit(full_stomachs$content_length.cm)))),
        cbind(full_stomachs$fork_length.cm,
              rep("predator", times = length(full_stomachs$fork_length.cm)))))

colnames(lengths) <- c("length", "hake")
lengths$length <- as.numeric(lengths$length)
  
# Histogram of lengths
length_hist <- ggplot(lengths, aes(x=length, fill=hake)) +
  geom_histogram() +
  stat_bin(binwidth = 3) +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("length (cm)") + ylab(" ")
length_hist

ggsave(filename="plots/diet/old/length_hist.png", length_hist,
       width=150, height=100, units="mm", dpi=300)

# Combine weights into one dataframe, excluding duplicate contents weights from
# multiple length observations
weights <- as.data.frame(
  rbind(cbind(unique((full_stomachs$content_weight.g) / 1000),  # convert to kg
              rep("prey", times = length(unique(full_stomachs$content_weight.g)))),
        cbind(full_stomachs$organism_weight.kg,  
              rep("predator", times = length(full_stomachs$organism_weight.kg)))))

colnames(weights) <- c("weight", "hake")
weights$weight <- as.numeric(weights$weight)

# Histogram of weights
weight_hist <- ggplot(weights, aes(x=weight, fill=hake)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("weight (kg)") + ylab(" ")
weight_hist

ggsave(filename="plots/diet/old/weight_hist.png", weight_hist,
       width=150, height=100, units="mm", dpi=300)


# Comparison with historic cannibalism data -----------------------------------
diet_history <- read.csv("data/diet/old/From Isaac/hake_diet_history.csv")
history2 <- melt(diet_history, id.vars = c("stage", "source", "year"))

history_plot <- ggplot(history2, aes(x=source, y=value, fill=stage)) +
  geom_bar(stat = "identity") +
  theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  ylab(" ") +
  facet_wrap(~variable, ncol=1, scales = "free_y")
history_plot

ggsave(filename="plots/diet/old/diet_history.png", history_plot,
       width=300, height=250, units="mm", dpi=300)
