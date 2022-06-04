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

# Calculate highest predation by relative weight
all_prey_wt <- pred_prey %>%
  group_by(Predator_Com_Name) %>%
  filter(!is.na(Prey_Weight_g)) %>%
  summarize(total_wt = sum(Prey_Weight_g))

hake_prey_wt <- pred_prey %>%
  filter(Prey_Com_Name == "Pacific Hake") %>%
  group_by(Predator_Com_Name) %>%
  filter(!is.na(Prey_Weight_g)) %>%
  summarize(hake_wt = sum(Prey_Weight_g))

hake_wt <- merge(all_prey_wt, hake_prey_wt, by = "Predator_Com_Name") %>%
  arrange(-hake_wt)

# Combine and plot
highest <- rbind(data.frame(predator = hake_wt$Predator_Com_Name[1:15], 
                            highest = hake_wt$hake_wt[1:15], 
                            variable = rep("relative weight (g)", 15)),
                data.frame(predator = high_n$Predator_Com_Name[1:15],
                           highest = high_n$freq[1:15], 
                           variable = rep("relative occurrence (top 15)", 15)))

hake_pred_plot <- ggplot(highest, aes(x = reorder(predator, highest), y = highest)) +
  geom_bar(position = "dodge", stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
  theme_sleek() +
  xlab(" ") + ylab(" ") +
  facet_wrap(~ variable, scales = "free")
hake_pred_plot

ggsave(filename = "plots/diet/top_hake_predators.png", hake_pred_plot, 
       width=300, height=100, units="mm", dpi=300)
