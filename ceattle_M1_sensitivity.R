# Sensitivity to M1 (estimated/fixed/prior)
# devtools::install_github("grantdadams/Rceattle", ref = "dev")
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggview)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

### Load in models - Rdata objects --------------------------------------------
load("models/ms_estM1.Rdata")
load("models/ms_fixM1.Rdata")
load("models/ms_priorM1.Rdata")

# Compare fits
model_fits <- rbind(cbind(model = "MS est M1", ms_estM1$fit),
                    cbind(model = "MS fix M1", ms_fixM1$fit),
                    cbind(model = "MS prior M1", ms_priorM1$fit))
model_summary <- cbind(ms_estM1$summary,
                       ms_fixM1$summary[, 3],
                       ms_priorM1$summary[, 3])[, -(1:2)]
colnames(model_summary) <- c("MS est M1", "MS fix M1", "MS prior M1")

### Plot  -----------------------------------------------------------
start_yr <- ms_estM1$model$data_list$styr
end_yr <- 2022
years <- start_yr:end_yr

# Helper function for extracting -by-age data from CEATTLE
extract_byage <- function(result, name, type) {
  df <- as.data.frame(as.table(result))
  
  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(years)
  colnames(df) <- c("age", "year", type)
  
  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"
  
  return(df)
}

ceattle_biomass <- function(run, name) {
  ssb <- (c(run$quantities$biomassSSB) * 2)
  biom <- c(run$quantities$biomass)
  biom_sd <- run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]
  wide <- as.data.frame(cbind(years, ssb, biom))
  colnames(wide) <- c("year", "SSB", "Total Biomass")
  all_biom <- melt(wide, id.vars = "year")
  all_biom2 <- cbind(all_biom, 
                     error = c(run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")], 
                               run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]), 
                     model = rep(name, length(all_biom$year)))  
  colnames(all_biom2)[2:3] <- c("type", "value")
  return(all_biom2)
}

biomass <- rbind.data.frame(ceattle_biomass(ms_estM1$model, "estimated"),
                            ceattle_biomass(ms_fixM1$model, "fixed"),
                            ceattle_biomass(ms_priorM1$model, "prior"))

recruitment <- data.frame(year = years,
                          estimated = c(ms_estM1$model$quantities$R),
                          fixed = c(ms_fixM1$model$quantities$R),
                          prior = c(ms_priorM1$model$quantities$R))
recruitment <- melt(recruitment, id.vars = "year", variable.name = "model")
recruitment$error <- c(ms_estM1$model$sdrep$sd[which(names(ms_estM1$model$sdrep$value) == "R")],
                       ms_fixM1$model$sdrep$sd[which(names(ms_fixM1$model$sdrep$value) == "R")],
                       ms_priorM1$model$sdrep$sd[which(names(ms_priorM1$model$sdrep$value) == "R")])
recruitment <- cbind.data.frame(year = recruitment$year,
                                type = rep("Recruitment"),
                                value = recruitment$value,
                                error = recruitment$error,
                                model = recruitment$model)

all_popdy <- rbind(biomass, recruitment)
all_popdy$value <- all_popdy$value / 1000000
all_popdy$error <- all_popdy$error / 1000000
all_popdy$type <- factor(all_popdy$type, 
                         labels = c("SSB (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))

# Add bounds for error & set 0 as minimum for plotting
all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
all_popdy$min[all_popdy$min < 0] <- 0
all_popdy$max <- all_popdy$value + (2 * all_popdy$error)

popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
  geom_line() +
  geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
  geom_vline(xintercept = 2019, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  ylim(0, NA) +
  ylab(" ") +
  labs(color = "model") +
  facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
  theme(strip.background = element_blank(), strip.placement = "outside") 
popdy_plot
