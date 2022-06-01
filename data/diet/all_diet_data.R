# Script for subsetting the whole diet database for predator hake

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(rgeos)
library(purrr)

path <- "data/diet/Full dataset/v4/"

combine_diet <- function(pred_species, prey_species) {
  # Read in and combine diet data by species ----------------------------------
  
  # All collection info
  collection <- read.csv(paste0(path, "collection_information_v4.csv"))
  
  # Filter predator dataset for just hake and collections from 1980 onwards
  predator <- read.csv(paste0(path, "predator_information_v4.csv")) %>%
    filter(Predator_Com_Name == pred_species) %>%
    filter(Collection_ID > 58) %>%  # Years before 1980
    filter(!is.na(FL_cm))  # Remove observations without fork lengths
  
  # Filter collection info by hake predator collection IDs
  collection_pred <- read.csv(paste0(path, "collection_information_v4.csv")) %>%
    filter(Collection_ID %in% predator$Collection_ID)
  
  # Filter prey composition dataset by hake predator IDs
  prey_comp <- read.csv(paste0(path, "prey_composition_v4.csv")) %>%
    filter(Predator_ID %in% predator$Predator_ID)
  
  # Filter prey size dataset by hake-predated-prey IDs
  prey_size <-  read.csv(paste0(path, "prey_size_v4.csv")) %>%
    filter(Prey_Comp_ID %in% prey_comp$Prey_Comp_ID)
  
  # Combine collection info with predator info
  all_pred <- merge(collection, predator, all.y = TRUE)[, c("Month", 
                                                            "Year", 
                                                            "Latitude", 
                                                            "Longitude",
                                                            "Predator_ID", 
                                                            "FL_cm", 
                                                            "Predator_Weight_kg")]
  
  
  all_prey <- merge(prey_comp, prey_size, all.y = TRUE)[, c("Prey_Com_Name", 
                                                            "Predator_ID", 
                                                            "Prey_Weight_g",
                                                            "Prey_Maturity",
                                                            "Prey_Length1")]
  
  # Combine pred & prey together
  all_stomachs <- merge(all_pred, all_prey, all.y = TRUE)
  
  predated <- all_stomachs %>%
    filter(Prey_Com_Name == prey_species)
  
  # Plot general trends in the data ------------------------------------------
  # Overall number of collections
  yearly_n <- collection %>%
    group_by(Year) %>%
    summarize(n = n()) %>%
    filter(Year >= 1980)
  
  collections_plot <- ggplot(yearly_n, aes(x = Year, y = n)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme_sleek() +
    xlab(" ") + ylab("Number of collections")
  
  # Top prey items by occurrence and weight
  high_wt <- prey_comp %>% 
    group_by(Prey_Com_Name) %>% 
    summarize(highest = sum(Prey_Weight_g)) %>%
    slice_max(n = 10, order_by = highest)
  
  high_n <- prey_comp %>% 
    group_by(Prey_Com_Name) %>% 
    summarize(highest = n()) %>%
    slice_max(n = 10, order_by = highest)
  
  highest <- rbind(cbind(high_wt, variable = rep("weight (top 10)", 10)), 
                   cbind(high_n, variable = rep("occurrence (top 10)", 10)))
  
  prey_sp_plot <- ggplot(highest, aes(x = Prey_Com_Name, y = highest,
                                      fill = ifelse(Prey_Com_Name == prey_species, "highlighted", "normal"))) +
    geom_bar(position = "dodge", stat = "identity", show.legend = FALSE) +
    coord_flip() +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
    theme_sleek() +
    xlab(" ") + ylab(" ") +
    facet_wrap(~ variable, scales = "free")
  
  return(list(all_pred, all_prey, predated, collections_plot, prey_sp_plot))
  
}

hake_hake <- combine_diet("Pacific Hake", "Pacific Hake")

# Look at plots
hake_hake[[4]]  # number of collections
hake_hake[[5]]  # top prey species

ggsave(filename = "plots/diet/hake_prey_species.png", hake_hake[[5]], 
       width=200, height=80, units="mm", dpi=300)


# Look at instances of cannibalism in hake ------------------------------------
# Label predator dataframe by instances of cannibalism & change order for plots
pred_type <- hake_hake[[1]]
pred_type$type <- ifelse(pred_type$Predator_ID %in% hake_hake[[3]]$Predator_ID, "cannibalistic", "general predator")
pred_type$type <- relevel(factor(pred_type$type), "general predator")

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


### Plot timing of sample collection ------------------------------------------
timing_all <- pred_type %>%
  group_by(Year, Month, type) %>%
  summarize(n = n()) %>%
  filter(!is.na(Year))

timing_yearly <- ggplot(timing_all, aes(x = as.factor(Month), y = n, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(limits=factor(1:12)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab("sampling month") + ylab(" ") +
  facet_wrap(~ Year)
timing_yearly

ggsave(filename = "plots/diet/timing_yearly.png", timing_yearly, 
       width=200, height=150, units="mm", dpi=300)

timing_overall <- ggplot(timing_all, aes(x = as.factor(Month), y = n, color = type, fill = type)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(limits=factor(1:12)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab("sampling month") + ylab(" ") 
timing_overall

ggsave(filename = "plots/diet/timing_overall.png", timing_overall, 
       width=170, height=100, units="mm", dpi=300)


### Plot location of sample collection ----------------------------------------
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
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) +  
  scale_x_continuous(breaks = seq(-135, -120, by = 10)) +
  scale_y_continuous(breaks = seq(35, 55, by = 10)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab(" ") + ylab(" ") +
  facet_wrap(~Year, ncol = 8)

ggsave(filename = "plots/diet/locations_yearly.png", location_yearly, 
       width=400, height=200, units="mm", dpi=300)

location_overall <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = location_all, aes(x = Longitude, y = Latitude, colour = type, size = n)) +
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) +
  scale_x_continuous(breaks = seq(-135, -115, by = 5)) +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab(" ") + ylab(" ") 
location_overall

ggsave(filename = "plots/diet/locations_overall.png", location_overall, 
       width=100, height=100, units="mm", dpi=300)


### Inset timing plots in yearly location plots -------------------------------
# Tutorial here: https://www.blopig.com/blog/2019/08/combining-inset-plots-with-facets-using-ggplot2/
get_inset <- function(df) {
  # Create plot for the inset 
  plot <- ggplot(df, aes(x = as.factor(Month), y = n, fill = type)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_x_discrete(limits = factor(1:12), breaks = c(1, 6, 12)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_sleek() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=rel(0.8)),  # inset axis tick font size
          plot.background = element_rect(fill='transparent', color=NA)) + # transparent so no overlap w/map
    theme(legend.position="none") 
  return(plot)
}

# Function for defining how the inset will be positioned
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

inset_plot <- get_inset(timing_all)  # Actually create insets

# How the insets will be mapped on to the main plots (applying above function)
insets <- timing_all %>% 
  split(f = .$Year) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.)), 
    data = data.frame(Year=unique(.$Year)),
    ymin = 30.7, ymax = 40, xmin = -141, xmax = -124))  # position of insets

# Bring everything together - add insets on to main plot (locations, created above)
location_timing <- location_yearly +
  coord_sf(xlim = c(-140, -115), ylim = c(31, 56), expand = FALSE) + 
  scale_x_continuous(breaks = seq(-135, -120, by = 10)) +
  insets
  
ggsave(filename = "plots/diet/location_timing.png", location_timing, 
       width=400, height=210, units="mm", dpi=300)
  
  
### Write predator & prey datasets to .csvs -----------------------------------
write.csv(hake_hake[[1]], "data/diet/Full dataset/full_hake_pred.csv")
write.csv(hake_hake[[2]], "data/diet/Full dataset/full_prey.csv")
