# Script for analyzing the hake diet data from US vessels. Data are split
# between the predators (all hake) and the prey that are hake.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(fishmethods)

stomach_summary <- read.csv("data/diet/hake_stomachs_summary.csv")
full_stomachs <- read.csv("data/diet/full_stomachs.csv")
hake_predator <- read.csv("data/diet/hake_predator_all.csv")
all_stomachs <- read.csv("data/diet/hake_stomachs_all.csv")

# Total stomachs per year -----------------------------------------------------
per_year <- ggplot(stomach_summary, aes(x=year)) +
  geom_histogram() +
  theme_sleek() +
  xlab("year") + ylab(" ")
per_year

ggsave(filename="plots/diet/stomachs_per_year.png", per_year,
       width=150, height=100, units="mm", dpi=300)

# Empty vs. hake_containing stomachs ------------------------------------------
empty_hake_length <- ggplot(all_stomachs, aes(x=fork_length.cm, fill=scientific_name)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, labels = c("empty", "hake")) +
  xlab("fork length (cm)") + ylab(" ") +
  labs(fill = "contents")
empty_hake_length

ggsave(filename="plots/diet/emptyvshake_length.png", empty_hake_length,
       width=150, height=100, units="mm", dpi=300)

empty_hake_age <- ggplot(all_stomachs, aes(x=age, fill=scientific_name)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, labels = c("empty", "hake")) +
  xlab("predator age") + ylab(" ") +
  labs(fill = "contents")
empty_hake_age

ggsave(filename="plots/diet/emptyvshake_age.png", empty_hake_age,
       width=150, height=100, units="mm", dpi=300)
  

# Predator hake dynamics ------------------------------------------------------
# Overall lengths 
pred_lengths <- ggplot(hake_predator, aes(x=fork_length.cm, fill=sex)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("fork length (cm)") + ylab(" ") 
pred_lengths

ggsave(filename="plots/diet/all_pred_lengths.png", pred_lengths,
       width=150, height=100, units="mm", dpi=300)

pred_age <- ggplot(hake_predator, aes(x=age, fill=sex)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("age") + ylab(" ")
pred_age

ggsave(filename="plots/diet/pred_age.png", pred_age,
       width=150, height=100, units="mm", dpi=300)

pred_sex <- ggplot(hake_predator, aes(x=year, fill=sex)) +
  geom_bar(stat="count") +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("year") + ylab(" ")
pred_sex

ggsave(filename="plots/diet/pred_sex.png", pred_sex,
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
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("weight (kg)") + ylab(" ")
weight_hist

ggsave(filename="plots/diet/weight_hist.png", weight_hist,
       width=150, height=100, units="mm", dpi=300)


# Comparison with historic cannibalism data -----------------------------------
diet_history <- read.csv("data/diet/From Isaac/hake_diet_history.csv")

history_plot <- ggplot(diet_history, aes(x=reorder(source, prop), y=prop, fill=stage)) +
  geom_bar(stat = "identity") +
  theme_sleek() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) 
history_plot

ggsave(filename="plots/diet/diet_history.png", history_plot,
       width=150, height=100, units="mm", dpi=300)




