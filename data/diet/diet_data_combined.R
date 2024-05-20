# Script for combining the SWFSC/CCTD data and the FEAT data from Alicia 
# together for Pacific Hake

library(dplyr)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(geos)
library(purrr)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

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

# Switch to age 20 accumulator age
FEAT_all <- FEAT_all %>% 
  filter(predator_age != "(blank)") %>%
  mutate(predator_age = as.numeric(predator_age)) %>%
  mutate(predator_age = ifelse(predator_age > 20, 20, predator_age))

FEAT_sampling <- FEAT_all %>%
  group_by(year, prey_name) %>%
  summarize(n = n()) %>%
  ggplot(aes(x = year, y = n, fill = prey_name)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  labs(fill = "prey species (n)") + ylab(" ")
FEAT_sampling

# ggsave(filename = "plots/diet/FEAT_sampling.png", FEAT_sampling, 
#        bg = "transparent", width=180, height=90, units="mm", dpi=300)


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

# Remove all prey ages for prey that isn't hake & remove 1980 (year w/ no ages)
all_data$prey_age[all_data$prey_name == "other"] <- NA

hake_data <- all_data %>%
  filter(prey_name == "Pacific Hake") 

# Look at instances of no age - all weights within year 1 range except for 1
no_age <- hake_data %>% filter(is.na(prey_length))
# Replace NAs with age 1
all_data$prey_age[is.na(all_data$prey_age) & all_data$prey_name == "Pacific Hake"] <- 1 

# Add estimate of correct age for very big hake prey
all_data$prey_age[all_data$prey_wt > 300 & all_data$prey_name == "Pacific Hake"] <- 5

# Rewmove 1980 w/ no ages
all_data <- all_data %>% filter(year != 1980)

# Write combined dataset
# write.csv(all_data, "data/diet/CCTD_FEAT_combined.csv", row.names = FALSE)


### Plot data -----------------------------------------------------------------
# Plot timing & location of instances of predation of interest
predation_all <- all_data %>%
  group_by(year, prey_name) %>%
  summarize(n = n()) %>%
  filter(!is.na(year))

predation_yearly <- ggplot(predation_all, aes(x = year, y = n, fill = prey_name)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  labs(fill = "prey species") + ylab("Stomachs (n)") + xlab("Year")
predation_yearly

# ggsave(filename = "plots/diet/hake_cannibalism.png", predation_yearly,
#        bg = "transparent", width=190, height=100, units="mm", dpi=300)

### Plot timing of sample collection ----------------------------------------
time_n <- all_data %>%
  group_by(year, month) %>%
  summarize(n_all = n()) %>%
  filter(!is.na(year))

time_n_cannibalism <- all_data %>%
  filter(prey_name == "Pacific Hake") %>%
  group_by(year, month) %>%
  summarize(n_cannibalism = n()) %>%
  filter(!is.na(year)) 

time_n_all <- left_join(time_n, time_n_cannibalism)
time_n_all$n_cannibalism[is.na(time_n_all$n_cannibalism)] <- 0
time_n_all$prop <- time_n_all$n_cannibalism / time_n_all$n_all
time_n_all <- time_n_all %>% arrange(year)
# time_n_all$month <- factor(time_n_all$month)

data_years <- c(1988:1992, 1995:1999, 2002, 2005, 2007, 2009, 2011:2013, 2015, 2017, 2019)
all_months <- cbind.data.frame(year = rep(data_years, each = 12),
                               month = rep(1:12, times = length(data_years)))
time_all_months <- left_join(all_months, time_n_all)
time_all_months[is.na(time_all_months)] <- 0

timing_yearly <- ggplot(time_all_months, aes(x = month, y = n_all, fill = prop)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = factor(1:12), breaks = c(2, 4, 6, 8, 10, 12)) +
  scale_fill_viridis(limits = c(0, 1), begin = 0.1, end = 0.9) +
  xlab("Sampling Month") + ylab(" ") +
  labs(fill = "Cannibalism Rate") +
  facet_wrap(~ year, ncol = 5)
timing_yearly

# ggsave(filename = "plots/diet/yearly_timing.png", timing_yearly,
#        bg = "transparent", width=200, height=140, units="mm", dpi=300)


### Plot location of sample collection --------------------------------------
loc_n <- all_data %>%
  group_by(year, latitude, longitude) %>%
  summarize(n_all = n()) %>%
  filter(!is.na(year))

loc_n_cannibalism <- all_data %>%
  filter(prey_name == "Pacific Hake") %>%
  group_by(year, latitude, longitude) %>%
  summarize(n_cannibalism = n()) %>%
  filter(!is.na(year))

loc_n_all <- left_join(loc_n, loc_n_cannibalism)
loc_n_all$n_cannibalism[is.na(loc_n_all$n_cannibalism)] <- 0
loc_n_all$prop <- loc_n_all$n_cannibalism / loc_n_all$n_all
loc_n_all <- loc_n_all %>% arrange(prop)

# New plot of cannibalism rate by year
annual_rate <- loc_n_all %>%
  group_by(year) %>%
  summarize(n = sum(n_all),
            prop = mean(prop)) %>%
  ggplot(., aes(x = year, y = n, fill = prop)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(limits = c(0, 0.4), begin = 0.1, end = 0.9) +
  ylab("Stomachs (n)") +
  labs(fill = "Proportion \nContaining \nHake") +
  theme(legend.title = element_text(hjust = 0))
annual_rate

ggsave(filename = "plots/diet/hake_cannibalism_rate.png", annual_rate,
       bg = "transparent", width=170, height=90, units="mm", dpi=300)

# Create a plot of location of observations by latitude and longitude
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry

locations <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = loc_n_all, aes(x = longitude, y = latitude, color = prop, size = n_all)) +
  coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) +
  scale_x_continuous(breaks = seq(-135, -115, by = 5)) +
  scale_y_continuous(breaks = seq(35, 55, by = 5)) +
  scale_color_viridis() +
  xlab(" ") + ylab(" ") + labs(color = "Proportion \nContaining \nHake", size = "Stomachs (n)") +
  facet_wrap(~year, ncol = 5)

### Inset timing plots in yearly location plots -----------------------------
# Tutorial here: https://www.blopig.com/blog/2019/08/combining-inset-plots-with-facets-using-ggplot2/
get_inset <- function(df) {
  # Create plot for the inset 
  plot <- ggplot(df, aes(x = factor(month), y = n_all, fill = prop)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(limits = factor(1:12), breaks = c(1, 12), labels = c("Jan", "Dec")) +
    scale_fill_viridis(limits = c(0, 1)) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=rel(0.9)),  # inset axis tick font size
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

inset_plot <- get_inset(time_n_all)  # Actually create insets

# How the insets will be mapped on to the main plots (applying above function)
insets <- time_all_months %>% 
  split(f = .$year) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.)), 
    data = data.frame(year=unique(.$year)),
    ymin = 30.7, ymax = 40, xmin = -141, xmax = -124))  # position of insets

# Bring everything together - add insets on to main plot (locations, created above)
location_timing <- locations +
  coord_sf(xlim = c(-140, -115), ylim = c(31, 56), expand = FALSE) + 
  scale_x_continuous(breaks = c(-135, -120)) +
  scale_y_continuous(breaks = c(32, 55)) +
  insets
location_timing

ggsave(filename = "~/Desktop/locations_timing.tiff", location_timing,
       width=170, height=180, units="mm", dpi=500)
