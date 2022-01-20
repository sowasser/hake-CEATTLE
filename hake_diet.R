# Script for analyzing the hake diet data from US vessels. Data are split
# between the predators (all hake) and the prey that are hake.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(fishmethods)

specimens <- read.csv("data/diet/hake_specimens.csv")
contents <- read.csv("data/diet/hake_contents.csv")
combined <- read.csv("data/diet/diet_combined.csv")

# Length & weight -------------------------------------------------------------
# Contents lengths 
contents_lengths <- ggplot(contents, aes(x=length)) +
  geom_histogram() +
  stat_bin(bins=4) +
  theme_sleek() +
  xlab("content length (cm)") + ylab(" ")

ggsave(filename="plots/diet/contents_length.png", contents_lengths,
       width=150, height=100, units="mm", dpi=300)

# Contents weights
contents_weights <- ggplot(contents, aes(x=content_weight)) +
  geom_histogram() +
  stat_bin(bins=5) +
  theme_sleek() +
  xlab("weight (g)") + ylab(" ")

ggsave(filename="plots/diet/contents_weight.png", contents_weights,
       width=150, height=100, units="mm", dpi=300)

# Specimen lengths
specimen_lengths <- ggplot(combined, aes(x=fork_length)) +
  geom_histogram() +
  stat_bin(bins=5) +
  theme_sleek() +
  xlab("specimen length (cm)") + ylab(" ")

ggsave(filename="plots/diet/specimen_lengths.png", specimen_lengths,
       width=150, height=100, units="mm", dpi=300)

# Specimen weights
specimen_weights <- ggplot(combined, aes(x=organism_weight_kg)) +
  geom_histogram() +
  stat_bin(bins=5) +
  theme_sleek() +
  xlab("organism weight (kg)") + ylab(" ")

ggsave(filename="plots/diet/specimen_weights.png", specimen_weights,
       width=150, height=100, units="mm", dpi=300)

# Estimate ages ---------------------------------------------------------------
# Compare assessment length-at-age to specimen & contents lengths
hake_maturity_data <- read.csv("~/Desktop/Local/hake-assessment-master/data/hake-maturity-data.csv")
length_age_all <- na.omit(hake_maturity_data[, 11:12])

# Check von Bertalanffy for estimating ages in hake diet data -----------------
# Estimates for Linf/Sinf & K from 2011 stock assessment
von_bert <- growth(size = length_age_all$Length_cm, age = length_age_all$Age,
                   Sinf = 50, K = 0.45, t0 = 0)
