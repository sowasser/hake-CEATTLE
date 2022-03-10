# Script for analyzing the hake diet data from US vessels. Data are split
# between the predators (all hake) and the prey that are hake.

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(fishmethods)
library(tidyr)

stomach_summary <- read.csv("data/diet/hake_stomachs_summary.csv")
full_stomachs <- read.csv("data/diet/full_stomachs.csv")
hake_predator <- read.csv("data/diet/hake_predator_all.csv")
all_stomachs <- read.csv("data/diet/hake_stomachs_all.csv")

# Total stomachs per year -----------------------------------------------------
per_year <- ggplot(stomach_summary, aes(x=year, y=num_stomachs_tow, fill=type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  xlab("year") + ylab("total stomachs")
per_year

ggsave(filename="plots/diet/stomachs_per_year.png", per_year,
       width=150, height=100, units="mm", dpi=300)

# Empty vs. hake_containing stomachs ------------------------------------------
empty_hake_length <- ggplot(hake_predator, aes(x=fork_length.cm, fill=scientific_name)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, labels = c("empty", "hake"), begin = 0.1, end = 0.9) +
  xlab("fork length (cm)") + ylab(" ") +
  labs(fill = "contents")
empty_hake_length

ggsave(filename="plots/diet/emptyvshake_length.png", empty_hake_length,
       width=150, height=100, units="mm", dpi=300)

empty_hake_age <- ggplot(hake_predator, aes(x=age, fill=scientific_name)) +
  geom_histogram() +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE, labels = c("empty", "hake"), begin = 0.1, end = 0.9) +
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
history2 <- melt(diet_history, id.vars = c("stage", "source", "year"))

history_plot <- ggplot(history2, aes(x=source, y=value, fill=stage)) +
  geom_bar(stat = "identity") +
  theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  ylab(" ") +
  facet_wrap(~variable, ncol=1, scales = "free_y")
history_plot

ggsave(filename="plots/diet/diet_history.png", history_plot,
       width=300, height=250, units="mm", dpi=300)

# Full analysis of hake predation ---------------------------------------------
# diet_all <- read.table("data/diet/From Isaac/hake_diet_history_all.csv", sep = ",", header = TRUE)
# diet2 <- melt(diet_all, id.vars = "Prey")
# diet3 <- diet2 %>% drop_na(value)
# 
# history_plot_all <- ggplot(diet3, aes(x=variable, fill=Prey)) +
#   geom_bar(position = "fill") +
#   theme_sleek() +
#   theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   ylab("proportion")
# 
# ggsave(filename="~/Desktop/diet_history.png", history_plot_all,
#        width=1000, height=1000, units="mm", dpi=300)

# Generalized diet
diet_sum <- read.csv("data/diet/From Isaac/hake_diet_history_sum.csv")
diet2 <- melt(diet_sum, id.vars = "Prey") %>% drop_na(value)

history_plot_sum <- ggplot(diet2, aes(x=variable, fill=Prey)) +
  geom_bar(position = "fill") +
  theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  ylab("proportion")

ggsave(filename="~/Desktop/diet_history_sum.png", history_plot_sum,
       width=700, height=500, units="mm", dpi=300)


