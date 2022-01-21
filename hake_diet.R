# Script for analyzing the hake diet data from US vessels. Data are split
# between the predators (all hake) and the prey that are hake.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(fishmethods)


specimens <- read.csv("data/diet/hake_specimens.csv")
contents <- read.csv("data/diet/hake_contents.csv")
combined <- read.csv("data/diet/diet_combined.csv")

# Length & weight -------------------------------------------------------------
lengths <- as.data.frame(
  rbind(cbind(na.omit(contents$length), 
              rep("contents", times = length(na.omit(contents$length)))),
        cbind(specimens$fork_length,
              rep("specimen", times = length(specimens$fork_length)))),
)

colnames(lengths) <- c("length", "hake")
lengths$length <- as.numeric(lengths$length)
  
# Histogram of lengths
length_hist <- ggplot(lengths, aes(x=length, fill=hake)) +
  geom_histogram() +
  stat_bin(binwidth = 3) +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75) +
  xlab("length (cm)") + ylab(" ")
# length_hist

ggsave(filename="plots/diet/length_hist.png", length_hist,
       width=150, height=100, units="mm", dpi=300)

# Combine weights into one dataframe, excluding duplicate contents weights from
# multiple length observations
weights <- as.data.frame(
  rbind(cbind(unique(contents$content_weight), 
              rep("contents", times = length(unique(contents$content_weight)))),
        cbind(specimens$organism_weight * 1000,  # convert kg to g
              rep("specimen", times = length(specimens$organism_weight)))),
)

colnames(weights) <- c("weight", "hake")
weights$weight <- as.numeric(weights$weight)

# Histogram of weights
weight_hist <- ggplot(weights, aes(x=weight, fill=hake)) +
  geom_histogram() +
  stat_bin(binwidth = 3) +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.25, end = 0.75) +
  xlab("weight (g)") + ylab(" ")
# weight_hist

ggsave(filename="plots/diet/weight_hist.png", weight_hist,
       width=150, height=100, units="mm", dpi=300)


# Estimate ages ---------------------------------------------------------------
hake_maturity_data <- read.csv("~/Desktop/Local/hake-assessment-master/data/hake-maturity-data.csv")
length_age_all <- na.omit(hake_maturity_data[, 11:12])

# Check von Bertalanffy for estimating ages in hake diet data
# Estimates for Linf/Sinf & K from 2011 stock assessment
von_bert <- growth(size = length_age_all$Length_cm, age = length_age_all$Age,
                   Sinf = 50, K = 0.45, t0 = 1)

# Apply inverted von Bertalanffy to length data 
inverted <- function(k, Linf, a0, length) {
  age <- c()
  for(l in length) {
    a <- -(1/k) * log(1 - (l / Linf)) + a0
    age <- c(age, a)
  }
  
  return(age)
}

estim_age_specimens <- inverted(k=0.45, Linf=50, a0=1, length=specimens$fork_length)
estim_age_contents <- inverted(k=0.45, Linf=50, a0=1, length=contents$length)
