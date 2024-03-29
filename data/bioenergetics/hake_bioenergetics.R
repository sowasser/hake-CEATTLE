# Script for investigating the sensitivity of the bioenergetics models to the
# parameters available for similar species - Walleye pollock and Atlantic cod.
# Data are from the Fish Bioenergetics 4.0 shiny app: http://fishbioenergetics.org/

library(reshape2)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

proxy_params <- read.csv("data/bioenergetics/proxy_bioen_params.csv")

eq_1 <- proxy_params[c(2:4), ]
eq_2 <- proxy_params[c(1, 5:6), ]
eq_2$Tco <- as.numeric(eq_2$Tco)
eq_2$Tcm <- as.numeric(eq_2$Tcm)

### Temperature-dependent consumption -----------------------------------------
temp_dependent <- function(Qc, Tco, Tcm, temp) {
  z <- (log(Qc) * (Tcm - Tco)) 
  y <- (log(Qc) * (Tcm - Tco + 2))
  
  x <- (((z^2) * (1 + (1 + (40/y))^0.5)^2) / 400)
  
  consumption <- c()
  for(i in temp) { 
    V <- ((Tcm - i) / (Tcm - Tco))
    rate <- ((V^x) * exp(x * (1 - V)))
    consumption <- c(consumption, rate)
  }
  
  return(consumption)
}

# Run equation with cod, adult, juvenile pollock parameters & hake estimates
temp_range <- seq(1, 20, by = 0.1)  # Hake thermal range from survey data
spp_temp_wide <- cbind(temp_dependent(eq_2[1, 5], eq_2[1, 6], eq_2[1, 7], temp_range), 
                       temp_dependent(eq_2[2, 5], eq_2[2, 6], eq_2[2, 7], temp_range), 
                       temp_dependent(2.5, 8, 14.5, temp_range),  # Hake estimates w/ temp from acoustic series
                       temp_dependent(2.5, 8, 10.5, temp_range))  # Hake estimates w/ kriged temp
colnames(spp_temp_wide) <- c("Atlantic cod - fb4", 
                             "Pollock - fb4", 
                             "Hake - survey temp",
                             "Hake - kriged temp")
spp_temp <- melt(as.data.frame(spp_temp_wide))
spp_temp <- cbind(spp_temp, temp = rep(temp_range, times=4))

# # Distinguish between literature values & estimated value for Hake (when that's ready)
# spp_temp <- cbind(spp_temp, ref = c(rep("a", times = (length(temp_range) * 3)), 
#                                     rep("b", times = length(temp_range))))
 
# Plot consumption rate                 
temp_rate <- ggplot(spp_temp, aes(x=temp, y=value)) +
  geom_line(aes(color=variable), linewidth=1) +
  # Following lines for distinguishing between lit & estimated hake values
  # geom_line(aes(color=variable, linetype=ref), size=1) +
  # scale_linetype_manual(values=c("longdash", "solid"), guide="none") +  
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.9) +  
  xlab("Temperature") + ylab("Specific Rate") +
  labs(color = "Species")
temp_rate

ggsave(filename="plots/bioenergetics/temp_consumption.png", temp_rate,
       bg = "transparent", width=180, height=90, units="mm", dpi=300)


# # Simpler version of the plot
# spp_temp_wide2 <- cbind(temp_dependent(eq_2[1, 5], eq_2[1, 6], eq_2[1, 7]), 
#                         temp_dependent(eq_2[2, 5], eq_2[2, 6], eq_2[2, 7]), 
#                         temp_dependent(2.5, 8, 10.5))  # Hake estimates w/ kriged temp
# colnames(spp_temp_wide2) <- c("Atlantic cod", 
#                               "Pollock", 
#                               "Pacific hake")
# spp_temp2 <- melt(as.data.frame(spp_temp_wide2))
# spp_temp2 <- cbind(spp_temp2, temp = rep(temp_range, times=3))
# temp_rate2 <- ggplot(spp_temp2, aes(x=temp, y=value)) +
#   geom_line(aes(color=variable), size=1) +
#   scale_color_viridis(discrete = TRUE, begin=0.1, end=0.9) +  
#   ylab("specific rate") +
#   labs(color = "species")
# temp_rate2
# 
# ggsave(filename="plots/bioenergetics/temp_consumption2.png", temp_rate2,
#        bg = "transparent", width=150, height=70, units="mm", dpi=300)



### Allometric mass function --------------------------------------------------
allometric_mass <- function(CA, CB, weights) {
  Cmax <- c()
  for(i in weights) {
    mass <- (CA * (i^CB))
    Cmax <- c(Cmax, mass)
  }
  return(Cmax)
}

# Run equation with cod, adult, juvenile pollock parameters & hake estimates
wt_range <- 10:400
spp_mass_wide <- cbind(allometric_mass(eq_2[1, 3], eq_2[1, 4], wt_range), 
                       allometric_mass(eq_2[2, 3], eq_2[2, 4], wt_range), 
                       allometric_mass(0.167, -0.460, wt_range),  # hake estimate from Francis (1983)
                       allometric_mass(0.0835, -0.460, wt_range))  # Francis (1983) estimate w/ CA/2  
colnames(spp_mass_wide) <- c("Atlantic cod - fb4", 
                             "Pollock (adult) - fb4", 
                             "Hake - Francis",
                             "Hake - Francis, CA/2")
spp_mass <- melt(as.data.frame(spp_mass_wide))
spp_mass <- cbind(spp_mass, weight = rep(wt_range, times=4))

# # Distinguish between literature values & estimated value for Hake (when that's ready)
# spp_mass <- cbind(spp_mass, ref = c(rep("a", times = (length(weights) * 5)), 
#                                     rep("b", times = length(weights))))

# Plot allometric mass
mass_rate <- ggplot(spp_mass, aes(x=weight, y=value)) +
  geom_line(aes(color=variable), size=1) +
  # Following lines for distinguishing between lit & estimated hake values
  # geom_line(aes(color=variable, linetype=ref), size=1) +
  # scale_linetype_manual(values=c("longdash", "solid"), guide="none") +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.9) + 
  ylab("Specific Rate") + xlab("Weight (g)") +
  labs(color = "Species")
mass_rate

ggsave(filename="plots/bioenergetics/allometric_mass.png", mass_rate,
       bg = "transparent", width=180, height=90, units="mm", dpi=300)


### Calculate relative foraging rate ------------------------------------------
# Calculate maximum consumption rate (CA * WCB * fT)
wtatage_kg <- read.csv("data/bioenergetics/wtatage_bioenergetics.csv")[, -1]  # mean annual weight-at-age
wtatage <- wtatage_kg * 1000
# Caclulate allometric mass function for annual weight-at-age
CA_wtatage <- t(apply(wtatage, 1, allometric_mass, CA=0.0835, CB=-0.460))  # by rows and transposed to match input

temp_yearly <- read.csv("data/temperature/ROMS_mean.csv")  # annual temperature
fT <- temp_dependent(2.5, 8, 10.5, temp_yearly$mean_temp)

multiply <- function(x) {
  out <- x * fT
}

max_rate <- apply(CA_wtatage, 2, multiply)
min(max_rate)
max(max_rate)

# Calculate observed consumption rate (for pollock)
pollock_obs <- as.data.frame(cbind(W = c(18.43, 49.09, 85.02, 100.53, 187.05, 210.77, 226.35, 254.81, 267.29, 302.36, 35.90),
                                   Cobs = c(0.022, 0.018, 0.015, 0.015, 0.008, 0.008, 0.007, 0.007, 0.008, 0.004, 0.027)))
ggplot(pollock_obs, aes(x = W, y = Cobs)) +
  geom_point()

# Estimate observed consumption rate for hake
summary(lm(Cobs ~ W, pollock_obs))
calc_Cobs <- function(wt) {
  out <- (-6.661e-05 * wt) + 2.316e-02
}

Cobs <- apply(wtatage, c(1, 2), calc_Cobs)

plot_cobs <- data.frame(cbind(Cobs = colMeans(Cobs), age = 1:15))
ggplot(plot_cobs, aes(x = age, y = Cobs)) +
  geom_point()

# Calculate relative foraging rate
RFR <- Cobs / max_rate

Pyrs <- RFR 
max(Pyrs)
min(Pyrs)


### Sensitivity of bioenergetics to temp anomaly ------------------------------
# Following the temperature anomalies determined from the hake survey data in
# this paper: https://www.int-res.com/abstracts/meps/v639/p185-197/
anomalies <- c(-1, -0.5, 0, 0.5, 1, 1.5, 2)
# Add mean temp from kriged temps where hake were present
anomaly_temp <- anomalies + 7.909423  

# Test sensitivity of temperature-dependent consumption curve
temp_dependent2 <- function(Qc, Tco, Tcm) {
  z <- (log(Qc) * (Tcm - Tco)) 
  y <- (log(Qc) * (Tcm - Tco + 2))
  
  x <- (((z^2) * (1 + (1 + (40/y))^0.5)^2) / 400)
  
  consumption <- c()
  for(i in anomaly_temp) { 
    V <- ((Tcm - i) / (Tcm - Tco))
    rate <- ((V^x) * exp(x * (1 - V)))
    consumption <- c(consumption, rate)
  }
  
  return(consumption)
}

temp_sen_wide <- cbind(temp_dependent2(2.6, 10, 15),  # Grant's pollock CEATTLE ex.
                       temp_dependent2(2.5, 8, 14.5),  # Hake estimates w/ temp from acoustic series
                       temp_dependent2(2.5, 8, 10.5))  # Hake estimates w/ kriged temp

colnames(temp_sen_wide) <- c("pollock - CEATTLE",
                             "hake estimate - survey temp",
                             "hake estimate - kriged temp")
temp_sen <- melt(as.data.frame(temp_sen_wide))
temp_sen <- cbind(temp_sen, anomaly = rep(anomalies, times=3))

# Plot resulting specific rates from the temp anomalies
temp_sensitivity <- ggplot(temp_sen, aes(x=anomaly, y=value)) +
  geom_point(aes(color=variable), size=2) +
  geom_line(aes(color=variable), size=1) +
  scale_color_viridis(discrete = TRUE, begin=0.1, end=0.9) + 
  ylab("specific rate") +
  xlab("temperature anomaly") +
  labs(color = "species")

ggsave(filename="plots/bioenergetics/temp_sensitivity.png", temp_sensitivity,
       bg = "transparent", width=200, height=100, units="mm", dpi=300)


# Plot difference between the 2 hake estimates for each anomaly
diff <- temp_dependent2(2.5, 8, 14.5) - temp_dependent2(2.5, 8, 10.5)
temp_diff <- data.frame(cbind(anomalies, diff))

rate_difference <- ggplot(temp_diff, aes(x=anomalies, y=diff)) +
  geom_point() +
  ylab("difference in specific rate") +
  xlab("temperature anomaly") 

ggsave(filename="plots/bioenergetics/rate_difference.png", rate_difference,
       bg = "transparent", width=200, height=100, units="mm", dpi=300)
