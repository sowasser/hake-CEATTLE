path <- "data/diet/Full dataset/v4/"

# All collection info
collection <- read.csv(paste0(path, "collection_information_v4.csv"))

predator <- read.csv(paste0(path, "predator_information_v4.csv"))

# Filter collection info by hake predator collection IDs
collection_pred <- read.csv(paste0(path, "collection_information_v4.csv")) %>%
  filter(Collection_ID %in% predator$Collection_ID)

# Filter prey composition dataset by hake predator IDs
prey_comp <- read.csv(paste0(path, "prey_composition_v4.csv"))

pred_prey <- merge(predator, prey_comp, by = "Predator_ID")

pred_prey_hake <- pred_prey %>%
  filter(Prey_Com_Name == "Pacific Hake") %>%

high_wt <- pred_prey %>%
  group_by(Predator_Com_Name, Prey_Com_Name) %>%
  summarize(total = sum(Prey_Weight_g)) %>%
  slice_max(n = 10, order_by = highest) %>%
  filter(Prey_Com_Name == "Pacific Hake")

high_n <- pred_prey %>%
  group_by(Predator_Com_Name, Prey_Com_Name) %>%
  summarize(highest = n()) %>%
  slice_max(n = 10, order_by = highest) %>%
  filter(Prey_Com_Name == "Pacific Hake") 

# Combine together and select *actual* top 10, in the case of a tie
highest <- rbind(cbind(high_wt[1:10,], variable = rep("weight (top 10)", 10)),
                 cbind(high_n[1:10,], variable = rep("occurrence (top 10)", 10)))

prey_sp_plot <- ggplot(highest, aes(x = reorder(Predator_Com_Name, highest), y = highest)) +
  geom_bar(position = "dodge", stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
  theme_sleek() +
  xlab(" ") + ylab(" ") +
  facet_wrap(~ variable, scales = "free")
prey_sp_plot
