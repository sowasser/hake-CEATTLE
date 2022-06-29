# Script for subsetting the whole diet database for predator hake

library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(rgeos)

path <- "data/diet/CCTD/v4/"

### Combine diet data, subset for species of interest -------------------------
combine_diet <- function(type, pred_species, prey_species, label_specific) {
  # All collection info
  collection <- read.csv(paste0(path, "collection_information_v4.csv"))
  
  # Filter predator dataset for just hake and collections from 1980 onwards
  if(type == "stomach") {
    predator <- read.csv(paste0(path, "predator_information_v4.csv")) %>%
      filter(Predator_Com_Name == pred_species) %>%
      filter(Collection_ID > 58) %>%  # Years before 1980
      filter(!is.na(FL_cm))  # Remove observations without fork lengths
  }
  if(type == "scat") {  # no lengths for scat data!
    predator <- read.csv(paste0(path, "predator_information_v4.csv")) %>%
      filter(Predator_Com_Name == pred_species) %>%
      filter(Collection_ID > 58) # Years before 1980
  }
  
  # Filter collection info by hake predator collection IDs
  collection_pred <- read.csv(paste0(path, "collection_information_v4.csv")) %>%
    filter(Collection_ID %in% predator$Collection_ID)
  
  # Filter prey composition dataset by hake predator IDs
  prey_comp <- read.csv(paste0(path, "prey_composition_v4.csv")) %>%
    filter(Predator_ID %in% predator$Predator_ID)
  
  if(pred_species == "Pacific Hake") {
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
  }
  if(pred_species == "Arrowtooth Flounder") {
    # No size info for arrowtooth, so combine predator info with prey comp
    all_pred <- merge(collection, predator, all.y = TRUE)[, c("Month", 
                                                              "Year", 
                                                              "Latitude", 
                                                              "Longitude",
                                                              "Predator_ID", 
                                                              "FL_cm", 
                                                              "Predator_Weight_kg")]
    all_prey <- prey_comp[, c("Prey_Com_Name", "Predator_ID")]
  }
  
  if(pred_species == "California Sea Lion") {
    # Filter prey size dataset by hake-predated-prey IDs
    # No predator size information for CA sea lion; only back-calculated length/weight for prey
    prey_size <-  read.csv(paste0(path, "prey_size_v4.csv")) %>%
      filter(Prey_Comp_ID %in% prey_comp$Prey_Comp_ID)
    
    # Combine collection info with predator info
    all_pred <- merge(collection, predator, all.y = TRUE)[, c("Month", 
                                                              "Year", 
                                                              "Latitude", 
                                                              "Longitude",
                                                              "Predator_ID")]
    all_prey <- merge(prey_comp, prey_size, all.y = TRUE)[, c("Prey_Com_Name", 
                                                              "Predator_ID", 
                                                              "Prey_Length_BC_mm",  
                                                              "Prey_Weight_Ind_BC_g")] 
  }
  
  # Combine pred & prey together
  all_stomachs <- merge(all_pred, all_prey, all.y = TRUE)
  
  predated <- all_stomachs %>%
    filter(Prey_Com_Name == prey_species)
  
  ### Plot general trends in the data -----------------------------------------
  # Top prey items by occurrence and weight
  high_wt <- prey_comp %>%
    group_by(Prey_Com_Name) %>%
    summarize(highest = sum(Prey_Weight_g)) %>%
    slice_max(n = 10, order_by = highest)

  high_n <- prey_comp %>%
    group_by(Prey_Com_Name) %>%
    summarize(highest = n()) %>%
    slice_max(n = 10, order_by = highest)

  # Combine together and select *actual* top 10, in the case of a tie
  highest <- rbind(cbind(high_wt[1:10,], variable = rep("weight (top 10)", 10)),
                   cbind(high_n[1:10,], variable = rep("occurrence (top 10)", 10)))

  prey_sp_plot <- ggplot(highest, aes(x = reorder(Prey_Com_Name, highest), y = highest,
                                      fill = ifelse(Prey_Com_Name == prey_species, "highlighted", "normal"))) +
    geom_bar(position = "dodge", stat = "identity", show.legend = FALSE) +
    coord_flip() +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, direction = -1) +
    theme_sleek() +
    xlab(" ") + ylab(" ") +
    facet_wrap(~ variable, scales = "free")
  
  # Plot timing & location of instances of predation of interest
  pred_type <- all_pred
  pred_type$type <- ifelse(pred_type$Predator_ID %in% predated$Predator_ID, label_specific, "general predator")
  pred_type$type <- relevel(factor(pred_type$type), "general predator")
  
  predation_all <- pred_type %>%
    group_by(Year, type) %>%
    summarize(n = n()) %>%
    filter(!is.na(Year))
  
  predation_yearly <- ggplot(predation_all, aes(x = Year, y = n, fill = type)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_sleek() +
    ylab(" ")
  
  ### Plot timing of sample collection ----------------------------------------
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
  
  ggsave(filename = "plots/diet/timing_yearly.png", timing_yearly, 
         width=200, height=150, units="mm", dpi=300)
  
  timing_overall <- ggplot(timing_all, aes(x = as.factor(Month), y = n, color = type, fill = type)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_x_discrete(limits=factor(1:12)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_sleek() +
    xlab("sampling month") + ylab(" ") 
  
  ggsave(filename = "plots/diet/timing_overall.png", timing_overall, 
         width=170, height=100, units="mm", dpi=300)
  
  ### Plot location of sample collection --------------------------------------
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
  
  location_overall <- ggplot(data = world) +
    geom_sf() +
    geom_point(data = location_all, aes(x = Longitude, y = Latitude, colour = type, size = n)) +
    coord_sf(xlim = c(-135, -115), ylim = c(31, 56), expand = FALSE) +
    scale_x_continuous(breaks = seq(-135, -115, by = 5)) +
    scale_y_continuous(breaks = seq(35, 55, by = 5)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_sleek() +
    xlab(" ") + ylab(" ") 
  
  ### Inset timing plots in yearly location plots -----------------------------
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
          geom = GeomCustomAnn,
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
      data = data.frame(Year = unique(.$Year)),
      ymin = 30.7, ymax = 40, xmin = -141, xmax = -124))  # position of insets
  
  # Bring everything together - add insets on to main plot (locations, created above)
  location_timing <- location_yearly +
    coord_sf(xlim = c(-140, -115), ylim = c(31, 56), expand = FALSE) + 
    scale_x_continuous(breaks = seq(-135, -120, by = 10)) +
    insets
  
  return(list(all_pred, all_prey, predated, pred_type,  # data
              prey_sp_plot, predation_yearly, location_overall, location_timing))  # plots
}


### Subset diet data for hake predator & prey ---------------------------------
hake_hake <- combine_diet(type = "stomach", "Pacific Hake", "Pacific Hake", "cannibalistic")

# Look at plots
hake_pred <- hake_hake[[1]]
hake_prey <- hake_hake[[2]]
hake_hake[[5]]  # top prey species
hake_hake[[6]]  # yearly predation by type
hake_hake[[7]]  # overall locations of predation by type

ggsave(filename = "plots/diet/hake_prey_species.png", hake_hake[[5]], 
       width=200, height=80, units="mm", dpi=300)
ggsave(filename = "plots/diet/hake_cannibalism.png", hake_hake[[6]], 
       width=170, height=100, units="mm", dpi=300)
ggsave(filename = "plots/diet/hake_locations_overall.png", hake_hake[[7]], 
       width=100, height=100, units="mm", dpi=300)
ggsave(filename = "plots/diet/hake_location_timing.png", hake_hake[[8]], 
       width=400, height=210, units="mm", dpi=300)


### Subset diet for arrowtooth flounder predator & hake prey ------------------
arrowtooth_hake <- combine_diet(type = "stomach", "Arrowtooth Flounder", "Pacific Hake", "hake predation")

ATF_pred <- arrowtooth_hake[[1]]
ATF_prey <- arrowtooth_hake[[2]]  
ATF_hake <- arrowtooth_hake[[3]]  # Almost no hake prey information for arrowtooth.

# Look at plots
arrowtooth_hake[[5]]  # top prey species
arrowtooth_hake[[6]]  # yearly predation by type
arrowtooth_hake[[7]]  # overall locations of predation by type

ggsave(filename = "plots/diet/ATF_prey_species.png", arrowtooth_hake[[5]], 
       width=200, height=80, units="mm", dpi=300)
ggsave(filename = "plots/diet/ATF_locations_overall.png", arrowtooth_hake[[7]], 
       width=100, height=100, units="mm", dpi=300)


### Subset of diet for CA sea lion predator & hake prey -----------------------
sealion_hake <- combine_diet(type = "scat", "California Sea Lion", "Pacific Hake", "hake predation")

CSL_pred <- sealion_hake[[1]]
CSL_prey <- sealion_hake[[2]]
CSL_hake <- sealion_hake[[3]]

# Look at plots
sealion_hake[[5]]  # top prey species
sealion_hake[[6]]  # yearly predation by type
sealion_hake[[7]]  # overall locations of predation by type

ggsave(filename = "plots/diet/CSL_prey_species.png", sealion_hake[[5]], 
       width=200, height=80, units="mm", dpi=300)
ggsave(filename = "plots/diet/CSL_locations_overall.png", sealion_hake[[7]], 
       width=100, height=100, units="mm", dpi=300)

CSL_hake_prey_size <- ggplot(CSL_hake, aes(x = Prey_Length_BC_mm / 10)) +
  ggdist::stat_slab(aes(thickness = stat(pdf*n)), scale = 0.7) +
  ggdist::stat_dotsinterval(side = "bottom", scale = 0.7, slab_size = NA) +
  theme_sleek() +
  xlab("prey hake length (cm)") + ylab(" ")
CSL_hake_prey_size

ggsave(filename = "plots/diet/CSL_hake_prey_size.png", CSL_hake_prey_size, 
       width=120, height=80, units="mm", dpi=300)

# Subset of diet for CA sea lion predator & ATF prey 
sealion_ATF <- combine_diet(type = "scat", "California Sea Lion", "Arrowtooth Flounder", "ATF predation")
CSL_prey_ATF <- sealion_ATF[[3]]  # no CA sea lion predation on ATF

  
### Write predator & prey datasets to .csvs -----------------------------------
write.csv(hake_hake[[1]], "data/diet/CCTD/hake_pred.csv", row.names = FALSE)
write.csv(hake_hake[[2]], "data/diet/CCTD/hake_prey.csv", row.names = FALSE)

write.csv(sealion_hake[[1]], "data/diet/CCTD/CSL_pred.csv", row.names = FALSE)
write.csv(sealion_hake[[2]], "data/diet/CCTD/CSL_prey.csv", row.names = FALSE)
