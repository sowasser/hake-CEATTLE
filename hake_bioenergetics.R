# Script for investigating the sensitivity of the bioenergetics models to the
# parameters available for similar species - Walleye pollock and Atlantic cod.
# Data are from the Fish Bioenergetics 4.0 shiny app: http://fishbioenergetics.org/

library(reshape2)
library(ggplot2)
library(ggsidekick)
library(viridis)

proxy_params <- read.csv("data/bioenergetics/proxy_bioen_params.csv")

eq_1 <- proxy_params[c(2:4), ]
eq_2 <- proxy_params[c(1, 5:6), ]
eq_2$Tco <- as.numeric(eq_2$Tco)
eq_2$Tcm <- as.numeric(eq_2$Tcm)

# Temperature-dependent consumption -------------------------------------------
temp_range <- seq(1, 20, by = 0.1)  # Hake thermal range from survey data

temp_dependent <- function(Qc, Tco, Tcm) {
  z <- (log(Qc) * (Tcm - Tco)) 
  y <- (log(Qc) * (Tcm - Tco + 2))
  
  x <- (((z^2) * (1 + (1 + (40/y))^0.5)^2) / 400)
  
  consumption <- c()
  for(i in temp_range) { 
    V <- ((Tcm - i) / (Tcm - Tco))
    rate <- ((V^x) * exp(x * (1 - V)))
    consumption <- c(consumption, rate)
  }
  
  return(consumption)
}

# Run equation with cod, adult, juvenile pollock parameters & hake estimates
spp_temp_wide <- cbind(temp_dependent(eq_2[1, 5], eq_2[1, 6], eq_2[1, 7]), 
                       temp_dependent(eq_2[2, 5], eq_2[2, 6], eq_2[2, 7]), 
                       temp_dependent(eq_2[3, 5], eq_2[3, 6], eq_2[3, 7]),
                       temp_dependent(2.6, 10, 15))  # Cq from Grant's pollock CEATTLE ex.
colnames(spp_temp_wide) <- c("Atlantic cod", 
                             "pollock (adult)", "pollock (juvenile)",
                             "hake - estimated")
spp_temp <- melt(as.data.frame(spp_temp_wide))
spp_temp <- cbind(spp_temp, temp = rep(temp_range, times=4), 
                  ref = c(rep("a", times = (length(temp_range) * 3)), 
                          rep("b", times = length(temp_range))))
                  
temp_rate <- ggplot(spp_temp, aes(x=temp, y=value)) +
  geom_line(aes(color=variable, linetype=ref), size=1) +
  scale_linetype_manual(values=c("longdash", "solid"), guide="none") +
  scale_color_viridis(discrete = TRUE) +  # invert colors
  theme_sleek() +
  ylab("specific rate") +
  labs(color = "species")

ggsave(filename="plots/bioenergetics/temp_consumption.pdf", temp_rate,
       width=200, height=100, units="mm", dpi=300)


# Allometric mass function ----------------------------------------------------
weights <- 10:400  # Range of hake weights
  
allometric_mass <- function(CA, CB) {
  Cmax <- c()
  for(i in weights) {
    mass <- (CA * (i^CB))
    Cmax <- c(Cmax, mass)
  }
  return(Cmax)
}

# Run equation with cod, adult, juvenile pollock parameters & hake estimates
spp_mass_wide <- cbind(allometric_mass(eq_2[1, 3], eq_2[1, 4]), 
                       allometric_mass(eq_2[2, 3], eq_2[2, 4]), 
                       allometric_mass(eq_2[3, 3], eq_2[3, 4]), 
                       allometric_mass(0.119009, -0.46024))  # values from Grant's pollock CEATTLE ex. 
colnames(spp_mass_wide) <- c("Atlantic cod", 
                             "pollock (adult)", "pollock (juvenile)",
                             "hake - estimated")
spp_mass <- melt(as.data.frame(spp_mass_wide))
spp_mass <- cbind(spp_mass, weight = rep(weights, times=4), 
                  ref = c(rep("a", times = (length(weights) * 3)), 
                          rep("b", times = length(weights))))

mass_rate <- ggplot(spp_mass, aes(x=weight, y=value)) +
  geom_line(aes(color=variable, linetype=ref), size=1) +
  scale_linetype_manual(values=c("longdash", "solid"), guide="none") +
  scale_color_viridis(discrete = TRUE) +  # invert colors
  theme_sleek() +
  ylab("specific rate") +
  labs(color = "species")

ggsave(filename="plots/bioenergetics/allometric_mass.pdf", mass_rate,
       width=200, height=100, units="mm", dpi=300)
