# Run CEATTLE for Pacific hake while testing different proportions of
# intraspecies predation, corrected with a Dirichlet multinomial distribution.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_221011.xlsx")

# # Run CEATTLE with the values as they are in the data file
# intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                  inits = NULL, # Initial parameters = 0
#                                  file = NULL, # Don't save
#                                  # debug = 1, # 1 = estimate, 0 = don't estimate
#                                  random_rec = FALSE, # No random recruitment
#                                  msmMode = 1, # Multispecies mode
#                                  phase = "default")


### Run CEATTLE with differing diet weight proportions ------------------------
# Read in different Dirichlet-corrected datasets
dirichlet_90s <- read.csv("data/diet/Dirichlet/Dirichlet_90s.csv")
dirichlet_recent <- read.csv("data/diet/Dirichlet/Dirichlet_recent.csv")


# Adapt weight proportions to replace those in the excel file & run CEATTLE
run_ceattle <- function(df, start, end, proj) {
  hake_intrasp$UobsWtAge <- df
  hake_intrasp$styr <- start
  hake_intrasp$endyr <- end
  hake_intrasp$projyr <- proj
  ceattle <- Rceattle::fit_mod(data_list = hake_intrasp,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # Multi-species mode
                               phase = "default")
  return(ceattle)
}

# Run low-cannibalism models
run_all <- run_ceattle(hake_intrasp$UobsWtAge, 1988, 2019, 2019)
run_90s <- run_ceattle(dirichlet_90s, 1988, 1999, 1999)
run_recent <- run_ceattle(dirichlet_recent, 2005, 2019, 2019)

# Run CEATTLE with no diet
nodiet_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                                inits = NULL, # Initial parameters = 0
                                file = NULL, # Don't save
                                # debug = 1, # 1 = estimate, 0 = don't estimate
                                random_rec = FALSE, # No random recruitment
                                msmMode = 0, # Single-species mode
                                phase = "default")


### Plot population dynamics --------------------------------------------------
ceattle_popdy <- function(run, name, years) {
  ssb <- (c(run$quantities$biomassSSB) * 2)
  biom <- c(run$quantities$biomass)
  wide <- as.data.frame(cbind(years, ssb, biom))
  colnames(wide) <- c("year", "SSB", "Total Biomass")
  all_biom <- melt(wide, id.vars = "year")
  all_biom2 <- cbind(all_biom, model = rep(name, length(all_biom$year)))  
  colnames(all_biom2)[2:3] <- c("type", "value")
  
  recruitment <- cbind(year = years,
                       type = rep("Recruitment"),
                       value = c(run$quantities$R),
                       model = name)
  
  popdy <- rbind(all_biom2, recruitment)
  
  return(popdy)
}

popdy <- rbind(ceattle_popdy(run_all, "all years", 1988:2019), 
               ceattle_popdy(run_90s, "1988-1999", 1988:1999), 
               ceattle_popdy(run_recent, "2005-2019", 2005:2019),
               ceattle_popdy(nodiet_run, "single-species", 1988:2019))

popdy$year <- as.numeric(popdy$year)
popdy$value <- as.numeric(popdy$value)
popdy$model <- factor(popdy$model,
                      levels = c("all years", "single-species", "1988-1999", "2005-2019"))

popdy_plot <- ggplot(popdy, aes(x=year, y=value, color = model, fill = model)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
  ylab(" ") +
  labs(color = "model") +
  facet_wrap(~type, ncol = 1, scales = "free_y")
popdy_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/dirichlet_popdy.png", popdy_plot, 
       bg = "transparent", width=170, height=140, units="mm", dpi=300)


### Plot mortality ------------------------------------------------------------
plot_mortality_custom <- function(Rceattle, file = NULL, incl_proj = FALSE, zlim = NULL, type = 0, width = 8,  height = 5.5, title = NULL, log = FALSE, minyr = NULL, theta = 155, species = NULL, maxage = NULL, title_cex = 10, M2 = TRUE) {
  
  # Convert single one into a list
  if(class(Rceattle) == "Rceattle"){
    Rceattle <- list(Rceattle)
  }
  
  if(length(Rceattle) > 1){
    stop("Can only plot one model")
  }
  
  # Extract data objects
  if(is.null(minyr)){ minyr <- Rceattle[[1]]$data_list$styr}
  
  Years <- minyr:Rceattle[[1]]$data_list$endyr
  if(incl_proj){
    Years <- minyr:Rceattle[[1]]$data_list$projyr
  }
  nyrs_vec <- length(Years)
  nyrs <- max(nyrs_vec)
  maxyr <- max((sapply(Years, max)))
  
  nspp <- Rceattle[[1]]$data_list$nspp
  spnames <- Rceattle[[1]]$data_list$spnames
  estdynamics <- Rceattle[[1]]$data_list$estDynamics
  nages <- Rceattle[[1]]$data_list$nages
  
  if(!is.null(maxage)){
    nages <- sapply(nages, function(x) ifelse(x > maxage, maxage, x))
  }
  
  minage <- Rceattle[[1]]$data_list$minage
  nsex <- Rceattle[[1]]$data_list$nsex
  
  # Get M
  M_array <-
    array(NA, dim = c(nspp, 2, max(nages), nyrs, length(Rceattle)))
  M1_array <-
    array(NA, dim = c(nspp, 2, max(nages), length(Rceattle)))
  for (i in 1:length(Rceattle)) {
    M1_array[, , ,i] <- Rceattle[[i]]$quantities$M1[,,1:max(nages)]
    if(!M2){
      M_array[, , , ,i] <- Rceattle[[i]]$quantities$M[,,1:max(nages),(1:nyrs)+(minyr - Rceattle[[1]]$data_list$styr)]
    }
    if(M2){
      M_array[, , , ,i] <- Rceattle[[i]]$quantities$M2[,,1:max(nages),(1:nyrs)+(minyr - Rceattle[[1]]$data_list$styr)]
    }
  }
  
  if(log){
    M1_array = log(M1_array)
    M_array = log(M_array)
  }
  
  # Plot limits
  zmax <- c()
  zmin <- c()
  for (i in 1:dim(M_array)[1]) {
    zmax[i] <- max(c(M_array[i,,,,], 0), na.rm = T)
    zmin[i] <- min(c(M_array[i,,,,], 0), na.rm = T)
  }
  
  # Plot trajectory
  loops <- ifelse(is.null(file), 1, 2)
  
  #################################
  # Mortality time series
  #################################
  if(is.null(species)){
    species <- 1:nspp
  }
  spp <- species
  
  # Species
  for(j in 1:nspp){
    sp <- j
    
    if(estdynamics[j] == 0 & sp %in% spp){
      
      # Sexes
      for(sex in 1:nsex[sp]){
        
        # Get sex for legend
        legend_sex = sex
        legend_sex2 = ifelse(sex == 1, "Female", "Male")
        if(nsex[sp] == 1){
          legend_sex <- 0
          legend_sex2 = "Combined"
        }
        
        # Save
        for (i in 1:loops) {
          if (i == 2) {
            filename <- paste0(file, "predation_and_residual_mortality_spp_",sp,"_sex_",legend_sex2,".png")
            png(
              file = filename ,
              width = width,
              height = height,
              units = "in",
              res = 300
            )
          }
          
          # Subset mortality data
          m_subset <- (M_array[j, sex, (1:nages[sp]), 1:nyrs, 1])
          
          # Get ages
          ages <- (1:(nages[sp])) - 1 + minage[sp]
          
          # Rearrange data
          data <- data.frame(Year = rep(Years, each = length(ages)), Age = rep(ages, length(Years)), M = c(m_subset))
          
          # Plot limits
          if(is.null(zlim)){
            zlim <- c(zmin[sp], zmax[sp])
          }
          
          if(is.character(zlim)){
            zlim <- c(min(zmin), max(zmax))
          }
          
          # Plot as tiles
          if(type == 0){
            p = ggplot2::ggplot(data, aes(y = Age, x = Year, zmin = zlim[1], zmax = zlim[2])) + 
              geom_tile(aes(fill = M))  + 
              scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) + 
              coord_equal() +  scale_x_continuous(expand = c(0, 0))+ 
              theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
            if(!is.null(title)){
              p = p + ggtitle(paste0(title,": ",spnames[j] )) + 
                theme(plot.title = element_text(size = title_cex))
            }
            
            scaleFUN <- function(x) sprintf("%.2f", x)  # set scaling function for legend
            
            if(log){
              p = p + scale_fill_viridis_c("log(M1 + M2)", limits = c(zlim[1], zlim[2]), labels = scaleFUN)
            } else {
              p = p + scale_fill_viridis_c("M1 + M2", limits = c(zlim[1], zlim[2]), labels = scaleFUN)
            }
            return(p)
          }
          
          # Plot as contours
          if(type == 1){
            print(ggplot2::ggplot(data, aes(y = Age, x = Year, z = M, zmin = zlim[1], zmax = zlim[2])) + geom_contour(colour = 1, size = 0.5) + geom_contour_filled()  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) +  scale_x_continuous(expand = c(0, 0)) + theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_fill_viridis_d("M1 + M2"))
          }
          
          # Plot as facets
          if(type == 2){
            p = ggplot(data=data, aes(x=Year, y = M, colour = Age, group = Age)) + theme( panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + geom_line(size = 2) + scale_color_viridis_c("Age")
            print(p)
          }
          
          # Plot as persp
          if(type == 3){
            par( mar=c(1 , 2 , 1 , 1) , tcl=-.25 , mgp=c(2 ,  1 ,  0) ,  oma=c(0 , 2 , 0 , 0))
            pmat = persp(y = Years, x = ages, z = m_subset, zlab = NA, zlim = zlim, xlab = "Age", ylab = "Year", theta = theta, ticktype = "detailed")
            mtext(ifelse(M2, "M2", "M"), side = 2, line = 0.5, at = 0)
            if(M2){
              text(-0.25,.15, labels = paste0("M1 = ",round((M1_array[j, sex, 1, 1]), 3)))
            }
            
            if(nsex[sp] == 1){
              mtext(paste0(title,": ",spnames[j]), side = 3, line = -2, at = 0)
            }
            if(nsex[sp] == 2){
              mtext(paste0(title,": ",spnames[j], " ",legend_sex2), side = 3, line = -2, at = 0)
            }
          }
          
          if (i == 2) {
            dev.off()
          }
        }
      }
    }
  }
}

# Combine mortality plots together & save
m_dirichlet <- gridExtra::grid.arrange(plot_mortality_custom(Rceattle = run_all, type = 0, title = "all years", maxage = 15),
                                       plot_mortality_custom(Rceattle = run_90s, type = 0, title = "1988-1999", maxage = 15) +
                                         scale_x_continuous(expand = c(0, 0), breaks = c(1988, 1992, 1996)),
                                       plot_mortality_custom(Rceattle = run_recent, type = 0, title = "2005-2019", maxage = 15),
                                       ncol = 2, nrow = 2, 
                                       layout_matrix = rbind(c(1,1), c(2,3)))

ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/M_dirichlet.png", 
       m_dirichlet, width=180, height = 180, units = "mm", dpi=300)

