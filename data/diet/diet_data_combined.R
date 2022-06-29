# Script for combining the SWFSC/CCTD data and the FEAT data from Alicia 
# together for Pacific Hake

library(dplyr)
library(ggplot2)
library(viridis)
library(ggsidekick)
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(rgeos)
library(purrr)

### Read in & update CCTD data from SWFSC -------------------------------------
CCTD_pred <- read.csv("data/diet/CCTD/hake_aged_pred.csv")
CCTD_prey <- read.csv("data/diet/CCTD/hake_aged_prey.csv")

# Combine predators and prey together, generalize prey names to "other"
CCTD_all <- merge(CCTD_pred, CCTD_prey, by = "Predator_ID") %>%
  mutate(prey_name = ifelse(Prey_Com_Name == "Pacific Hake", "Pacific Hake", "other")) %>%
  select(Predator_ID, Month, Year, Latitude, Longitude, pred_ages, prey_name, 
         Prey_Weight_g, Prey_Length1, prey_ages) %>%
  filter(Year <= 2004)  # remove years covered by FEAT 

# Add column for source of data
CCTD_all <- cbind(CCTD_all, source = rep("CCTD", length(CCTD_all[, 1])))


### Read in & update FEAT data from NWFSC -------------------------------------
FEAT_data <- read.csv("data/diet/FEAT/all_FEAT_diet.csv")
FEAT_prey_lengths <- read.csv("data/diet/FEAT/hake_prey_lengths_FEAT.csv")

# Add prey lengths based on stomach ID
FEAT_all <- merge(FEAT_data, FEAT_prey_lengths, by = "stomach_uuid", all.x = TRUE)

# Expand to month and year columns
FEAT_all$tow_timestamp <- as.Date(FEAT_all$tow_timestamp)
FEAT_all$month <- lubridate::month(lubridate::ymd(FEAT_all$tow_timestamp))
FEAT_all$year <- lubridate::year(lubridate::ymd(FEAT_all$tow_timestamp))

# Generalize prey categories to Pacific Hake & "other", select columns
FEAT_all <- FEAT_all %>%
  mutate(prey_name = ifelse(prey_category_long == "Gadiformes", "Pacific Hake", "other")) %>%
  select(stomach_uuid, month, year, tow_latitude, tow_longitude, predator_age, 
         prey_name, content_wt_g, measure_value)
  

### Combine datasets ----------------------------------------------------------
# Add column for prey ages. All are age 1 according to calculations in the 
# diet_aging.R script
FEAT_all <- cbind(FEAT_all, 
                  prey_ages = rep(1, length(FEAT_all[, 1])),
                  source = rep("FEAT", length(FEAT_all[, 1])))

labels <- c("Predator_ID", "month", "year", "latitude", "longitude", "predator_age", 
            "prey_name", "prey_wt", "prey_length", "prey_age", "source")

colnames(CCTD_all) <- labels
colnames(FEAT_all) <- labels

all_data <- rbind(CCTD_all, FEAT_all)

# Remove all prey ages for prey that isn't hake
all_data$prey_age[all_data$prey_name == "other"] <- NA

hake_data <- all_data %>%
  filter(prey_name == "Pacific Hake")

write.csv(all_data, "data/diet/CCTD_FEAT_combined.csv", row.names = FALSE)


### Plot data -----------------------------------------------------------------
# Plot timing & location of instances of predation of interest
predation_all <- all_data %>%
  group_by(year, prey_name) %>%
  summarize(n = n()) %>%
  filter(!is.na(year))

predation_yearly <- ggplot(predation_all, aes(x = year, y = n, fill = prey_name)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  labs(fill = "prey species (n)") + ylab(" ")
predation_yearly

### Plot timing of sample collection ----------------------------------------
timing_all <- all_data %>%
  group_by(year, month, prey_name) %>%
  summarize(n = n()) %>%
  filter(!is.na(year))

timing_yearly <- ggplot(timing_all, aes(x = as.factor(month), y = n, fill = prey_name)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(limits=factor(1:12)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab("sampling month") + ylab(" ") +
  facet_wrap(~ year)
timing_yearly

timing_overall <- ggplot(timing_all, aes(x = as.factor(month), y = n, color = prey_name, fill = prey_name)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(limits=factor(1:12)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab("sampling month") + ylab(" ") + 
  labs(fill = "prey species (n)", color = "prey species (n)")
timing_overall


### Plot location of sample collection --------------------------------------
location_all <- all_data %>%
  group_by(year, latitude, longitude, prey_name) %>%
  summarize(n = n()) %>%
  filter(!is.na(year)) %>%
  arrange(prey_name)

# Create a plot of location of observations by latitude and longitude
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry

location_yearly <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = location_all, aes(x = longitude, y = latitude, colour = prey_name, size = n)) +
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) +  
  scale_x_continuous(breaks = seq(-135, -120, by = 10)) +
  scale_y_continuous(breaks = seq(35, 55, by = 10)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab(" ") + ylab(" ") +
  facet_wrap(~year, ncol = 7)

location_overall <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = location_all, aes(x = longitude, y = latitude, colour = prey_name, size = n)) +
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) +
  scale_x_continuous(breaks = seq(-135, -115, by = 5)) +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_sleek() +
  xlab(" ") + ylab(" ") + labs(color = "prey species")
location_overall

### Inset timing plots in yearly location plots -----------------------------
# Tutorial here: https://www.blopig.com/blog/2019/08/combining-inset-plots-with-facets-using-ggplot2/
get_inset <- function(df) {
  # Create plot for the inset 
  plot <- ggplot(df, aes(x = as.factor(month), y = n, fill = prey_name)) +
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
  split(f = .$year) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.)), 
    data = data.frame(year=unique(.$year)),
    ymin = 30.7, ymax = 40, xmin = -141, xmax = -124))  # position of insets

# Bring everything together - add insets on to main plot (locations, created above)
location_timing <- location_yearly +
  coord_sf(xlim = c(-140, -115), ylim = c(31, 56), expand = FALSE) + 
  scale_x_continuous(breaks = seq(-135, -120, by = 10)) +
  insets

location_timing
