# Script for analyzing the hake diet data from US vessels. Data are split
# between the predators (all hake) and the prey that are hake.

library(dplyr)
library(ggplot2)
library(ggsidekick)

specimens <- read.csv("data/diet/hake_specimens.csv")
contents <- read.csv("data/diet/hake_contents.csv")


# Stomach contents only -------------------------------------------------------
# Prey lengths 
lengths <- ggplot(prey, aes(x=length)) +
  geom_histogram() +
  theme_sleek() +
  xlab("content length (cm)")

ggsave(filename="plots/diet/contents_length.png", lengths,
       width=150, height=100, units="mm", dpi=300)

# Prey weights
weights <- ggplot(prey, aes(x=content_weight)) +
  geom_histogram() +
  theme_sleek() +
  xlab("weight (g)")

ggsave(filename="plots/diet/contents_length.png", lengths,
       width=150, height=100, units="mm", dpi=300)
