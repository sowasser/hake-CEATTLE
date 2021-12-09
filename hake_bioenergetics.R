# Script for investigating the sensitivity of the bioenergetics models to the
# parameters available for similar species - Walleye pollock and Atlantic cod.
# Data are from the Fish Bioenergetics 4.0 shiny app: http://fishbioenergetics.org/

proxy_params <- read.csv("data/bioenergetics/proxy_bioen_params.csv")

# Temperature-dependent growth function 
temp_dependent <- function(Qc, Tco, Tcm, Temp){
  z <- (log(Qc) * (Tcm - Tco)) 
  y <- (log(Qc) * (Tcm - Tco + 2))
  
  x <- (((z^2) * (1 + (1 + (40/y))^0.5)^2) / 400)
  
  V <- ((Tcm - Temp) / (Tcm - Tco))
  
  growth <- ((V^x) * exp(x * (1 - V)))
    
  return(growth)
}

temp_dependent(2, 10, 15, 11)
