# Which predators eat the most hake?
library(dplyr)
library(ggplot2)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

path <- "data/diet/CCTD/v4/"

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

# Plot highest predators
hake_pred_plot <- ggplot(high_n, aes(x = reorder(Predator_Com_Name, freq), y = freq)) +
  geom_bar(position = "dodge", stat = "identity", show.legend = FALSE, fill = "#482577") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  xlab(" ") + ylab("Relative Frequency of Hake Predation") 
hake_pred_plot

ggsave(filename = "plots/diet/hake_predators.png", hake_pred_plot, 
       bg = "transparent", width=180, height=100, units="mm", dpi=300)

