# Script for analyzing the hake diet data from US vessels. Data are split
# between the predators (all hake) and the prey that are hake.

library(dplyr)

predator <- read.csv("data/diet/hake_predator_05-19.csv")
prey <- read.csv("data/diet/hake_prey_19.csv")