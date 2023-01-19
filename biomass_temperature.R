# Script for looking at the comparison between trends in biomass (with and
# without cannibalism) and annual mean temperature.

library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggridges)
library(ggview)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

years <- 1988:2019

# Read in and transform data --------------------------------------------------
ss3_ssb <- cbind(read.table("data/assessment/ssb.txt")[23:54, 2:3], type = rep("SSB"))
ss3_biomass <- cbind(read.table("data/assessment/biomass.txt")[23:54, 2:3], type = rep("Total Biomass"))

biomass <- read.csv("data/ceattle_intrasp_biomass.csv")
nodiet_biom <- read.csv("data/ceattle_nodiet_biomass.csv")

ratio <- as.data.frame(cbind(year = years,
                             assessment = ss3_ssb[, 1]/ss3_biomass[, 1],
                             cannibalism = biomass[1:length(years), 3] / biomass[(length(years)+1):(length(years)*2), 3],
                             single_species = nodiet_biom[1:length(years), 3] / nodiet_biom[(length(years)+1):(length(years)*2), 3]))

colnames(ratio) <- c("year", "Assessment", "CEATTLE - cannibalism", "CEATTLE - single-species")
ratio <- melt(ratio, id.vars = "year", variable.name = "model")
ratio$type <- "SSB/Biomass"

### Plot biomass and temperature together -------------------------------------
# Read in temp data from ROMS
roms <- cbind(melt(read.csv("data/temperature/ROMS_mean.csv"), id.vars = "Year"),
              type = rep("Mean Temperature")) %>%
  filter(Year >= 1988 & Year <= 2019) 
colnames(roms)[1:2] <- c("year", "model")
roms$model <- rep("ROMS")

bio_temp <- rbind(ratio, roms)
bio_temp$model <- factor(bio_temp$model, levels = c("ROMS", "Assessment", 
                                                    "CEATTLE - single-species", 
                                                    "CEATTLE - cannibalism"))

bio_temp_plot <- ggplot(bio_temp, aes(x=year, y=value)) +
  geom_line(aes(color=model)) +
  scale_color_manual(values = c("black", "#badb29", "#2b9089", "#4c2972")) +
  ylab(" ") +
  facet_wrap(~type, scales = "free_y", ncol = 1)
bio_temp_plot

ggsave(filename = "plots/CEATTLE/biomass_temp.png", bio_temp_plot,
       width = 150, height = 100, units = "mm")