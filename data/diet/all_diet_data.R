# Script for subsetting the whole diet database for predator hake

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(rgeos)

path <- "data/diet/Full dataset/v4/"

### Read in and subset datasets -----------------------------------------------
# Filter predator dataset for just hake and collections from 1980 onwards
predator_hake <- read.csv(paste0(path, "predator_information_v4.csv")) %>%
  filter(Predator_Com_Name == "Pacific Hake") %>%
  filter(Collection_ID > 58) %>%  # Years before 1980
  filter(!is.na(FL_cm))  # Remove observations without fork lengths

# Filter collection info by hake predator collection IDs
collection_hake <- read.csv(paste0(path, "collection_information_v4.csv")) %>%
  filter(Collection_ID %in% predator_hake$Collection_ID)

# Filter prey composition dataset by hake predator IDs
prey_of_hake_comp <- read.csv(paste0(path, "prey_composition_v4.csv")) %>%
  filter(Predator_ID %in% predator_hake$Predator_ID)

# Filter prey size dataset by hake-predated-prey IDs
prey_of_hake_size <-  read.csv(paste0(path, "prey_size_v4.csv")) %>%
  filter(Prey_Comp_ID %in% prey_of_hake_comp$Prey_Comp_ID)


### Combine then refine datasets ----------------------------------------------
# Combine collection info with predator info
all_pred <- merge(collection_hake, predator_hake, all.y = TRUE)[, c("Month", 
                                                                    "Year", 
                                                                    "Latitude", 
                                                                    "Longitude",
                                                                    "Predator_ID", 
                                                                    "FL_cm", 
                                                                    "Predator_Weight_kg")]

# Combine prey-of-hake datasets together
all_prey <- merge(prey_of_hake_comp, prey_of_hake_size, all.y = TRUE)[, c("Prey_Com_Name", 
                                                                          "Predator_ID", 
                                                                          "Prey_Weight_g",
                                                                          "Prey_Maturity",
                                                                          "Prey_Length1")]


# Combine pred & prey togehter - this is just full hake stomachs.
full_stomachs <- merge(all_pred, all_prey, all.y = TRUE)

# Cannibalism instances
hake_hake <- full_stomachs %>%
  filter(Prey_Com_Name == "Pacific Hake")

# Label predator dataframe by instances of cannibalism & change order for plots
pred_type <- all_pred
pred_type$type <- ifelse(pred_type$Predator_ID %in% hake_hake$Predator_ID, "cannibalistic", "general predator")
pred_type$type <- relevel(factor(pred_type$type), "general predator")


### Plot general trends in data -----------------------------------------------
# Top prey items by occurrence and weight
high_wt <- prey_of_hake_comp %>% 
  group_by(Prey_Com_Name) %>% 
  summarize(highest = sum(Prey_Weight_g)) %>%
  slice_max(n = 10, order_by = highest)

high_n <- prey_of_hake_comp %>% 
  group_by(Prey_Com_Name) %>% 
  summarize(highest = n()) %>%
  slice_max(n = 10, order_by = highest)

highest <- rbind(cbind(high_wt, variable = rep("weight (top 10)", 10)), 
                 cbind(high_n, variable = rep("occurrence (top 10)", 10)))

prey_sp <- ggplot(highest, aes(x = Prey_Com_Name, y = highest,
                               fill = ifelse(Prey_Com_Name == "Pacific Hake", "highlighted", "normal"))) +
  geom_bar(position = "dodge", stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
  theme_sleek() +
  xlab(" ") + ylab(" ") +
  facet_wrap(~ variable, scales = "free")
prey_sp

ggsave(filename = "plots/diet/hake_prey_species.png", prey_sp, 
       width=200, height=80, units="mm", dpi=300)

# Lengths
pred_fl <- ggplot(predator_hake, aes(x = FL_cm)) +
  geom_histogram() +
  xlab("fork length (cm)") + ylab(" ") +
  theme_sleek()
pred_fl

prey_length <- ggplot(hake_hake, aes(x = (Prey_Length1/10))) +
  geom_histogram() +
  xlab("length (cm)") + ylab(" ") +
  theme_sleek()
prey_length

# Predator hake stomachs including prey hake by year
cannibalism_all <- pred_type %>%
  group_by(Year, type) %>%
  summarize(n = n()) %>%
  filter(!is.na(Year))

cannibalism <- ggplot(cannibalism_all, aes(x = Year, y = n, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  ylab(" ")
cannibalism

ggsave(filename = "plots/diet/cannibalism.png", cannibalism, 
       width=170, height=100, units="mm", dpi=300)

# Timing
timing_all <- pred_type %>%
  group_by(Year, Month, type) %>%
  summarize(n = n()) %>%
  filter(!is.na(Year))

timing_yearly <- ggplot(timing_all, aes(x = as.factor(Month), y = n, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab("sampling month") + ylab(" ") +
  facet_wrap(~ Year)
timing_yearly

ggsave(filename = "plots/diet/timing_yearly.png", timing_yearly, 
       width=200, height=150, units="mm", dpi=300)


timing_overall <- ggplot(timing_all, aes(x = as.factor(Month), y = n, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab("sampling month") + ylab(" ") 
timing_overall

ggsave(filename = "plots/diet/timing_overall.png", timing_overall, 
       width=170, height=100, units="mm", dpi=300)


# Location
location_all <- pred_type %>%
  group_by(Year, Latitude, Longitude, type) %>%
  summarize(n = n()) %>%
  filter(!is.na(Year)) %>%
  arrange(type)

# Create a plot of location of observations by latitude and longitude
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry

location_yearly <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = location_all, aes(x = Longitude, y = Latitude, colour = type, size = n)) +
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) + xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(breaks = seq(-130, -110, by = 10)) +
  scale_y_continuous(breaks = seq(30, 50, by = 10)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab(" ") + ylab(" ") +
  facet_wrap(~Year, ncol = 7)

ggsave(filename = "plots/diet/locations_yearly.png", location_yearly, 
       width=300, height=300, units="mm", dpi=300)

location_overall <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = location_all, aes(x = Longitude, y = Latitude, colour = type, size = n)) +
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) + xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(breaks = seq(-135, -115, by = 5)) +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab(" ") + ylab(" ") 
location_overall

ggsave(filename = "plots/diet/locations_overall.png", location_overall, 
       width=100, height=100, units="mm", dpi=300)
  
  
### Write predator & prey datasets to .csvs -----------------------------------
write.csv(all_pred, "data/diet/Full dataset/full_hake_pred.csv")
write.csv(all_prey, "data/diet/Full dataset/full_prey.csv")
