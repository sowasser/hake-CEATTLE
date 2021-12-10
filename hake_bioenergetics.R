# Script for investigating the sensitivity of the bioenergetics models to the
# parameters available for similar species - Walleye pollock and Atlantic cod.
# Data are from the Fish Bioenergetics 4.0 shiny app: http://fishbioenergetics.org/

proxy_params <- read.csv("data/bioenergetics/proxy_bioen_params.csv")

# Temperature-dependent consumption -------------------------------------------
temp_range <- 5:15  # Hake thermal range from survey data

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

test <- temp_dependent(2, 10, 15)

# Allometric mass function ----------------------------------------------------
weights <- 100:400  # Range of hake weights
  
allometric_mass <- function(CA, CB) {
  Cmax <- c()
  for(i in weights) {
    mass <- (CA * (i^CB))
    Cmax <- c(Cmax, mass)
  }
  return(Cmax)
}

test2 <- allometric_mass(0.3, -0.3)
