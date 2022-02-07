# Script for analyzing the hake diet data from US vessels. Data are split
# between the predators (all hake) and the prey that are hake.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(fishmethods)

stomach_summary <- read.csv("data/diet/hake_stomachs_summary.csv")
full_stomachs <- read.csv("data/diet/full_stomachs.csv")

# Total stomachs per year -----------------------------------------------------
per_year <- ggplot(stomach_summary, aes(x=year)) +
  geom_histogram() +
  # stat_bin(binwidth = 3) +
  theme_sleek() +
  xlab("year") + ylab(" ")
per_year


# Length & weight -------------------------------------------------------------
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
  scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75) +
  xlab("length (cm)") + ylab(" ")
length_hist

ggsave(filename="plots/diet/length_hist.png", length_hist,
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
  scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.7) +
  xlab("weight (kg)") + ylab(" ")
weight_hist

ggsave(filename="plots/diet/weight_hist.png", weight_hist,
       width=150, height=100, units="mm", dpi=300)
