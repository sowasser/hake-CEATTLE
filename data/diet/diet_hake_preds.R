# Which predators eat the most hake?

library(dplyr)
library(ggplot2)
library(ggsidekick)

path <- "data/diet/Full dataset/v4/"

predator <- read.csv(paste0(path, "predator_information_v4.csv"))
prey_comp <- read.csv(paste0(path, "prey_composition_v4.csv"))

pred_prey <- merge(predator, prey_comp, by = "Predator_ID")

high_n <- pred_prey %>%
  group_by(Predator_Com_Name, Prey_Com_Name) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(Prey_Com_Name == "Pacific Hake") %>%
  arrange(-freq)

prey_sp_plot <- ggplot(high_n, aes(x = reorder(Predator_Com_Name, freq), y = freq)) +
  geom_bar(position = "dodge", stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
  theme_sleek() +
  xlab(" ") + ylab(" ")
prey_sp_plot
