# Script for analyzing the hake diet data from US vessels. Data are split
# between the predators (all hake) and the prey that are hake.

library(dplyr)
library(ggplot2)
library(ggsidekick)

predator <- read.csv("data/diet/hake_predator_05-19.csv")
prey <- read.csv("data/diet/hake_prey_19.csv")


# Prey ------------------------------------------------------------------------
# Prey lengths 
lengths <- ggplot(prey, aes(x=length)) +
  geom_histogram() +
  theme_sleek() +
  xlab("content weight (g)")

lengths

# Prey weights
weights <- ggplot(prey, aes(x=content_weight)) +
  geom_histogram() +
  theme_sleek() +
  xlab("length (cm)")
weights
