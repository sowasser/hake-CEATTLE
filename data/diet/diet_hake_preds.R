# Which predators eat the most hake?

library(dplyr)
library(ggplot2)
library(ggsidekick)

path <- "data/diet/Full dataset/v4/"

# Read in all predator and prey data
predator <- read.csv(paste0(path, "predator_information_v4.csv"))
prey_comp <- read.csv(paste0(path, "prey_composition_v4.csv"))

pred_prey <- merge(predator, prey_comp, by = "Predator_ID")

# Calculate highest predation by relative occurrence
high_n <- pred_prey %>%
  group_by(Predator_Com_Name, Prey_Com_Name) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(Prey_Com_Name == "Pacific Hake") %>%
  arrange(-freq)

# Combine and plot
highest <- data.frame(predator = high_n$Predator_Com_Name[1:15],
                      highest = high_n$freq[1:15])

hake_pred_plot <- ggplot(high_n, aes(x = reorder(Predator_Com_Name, freq), y = freq)) +
  geom_bar(position = "dodge", stat = "identity", show.legend = FALSE) +
  coord_flip() +
  theme_sleek() +
  xlab(" ") + ylab("Relative frequency of hake predation") 
hake_pred_plot

ggsave(filename = "plots/diet/hake_predators.png", hake_pred_plot, 
       width=300, height=100, units="mm", dpi=300)
