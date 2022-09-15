# Run CEATTLE for Pacific hake while testing different proportions of
# intraspecies predation

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())

hake_intrasp <- Rceattle::read_data( file = "data/hake_intrasp_220914.xlsx")

# # Run CEATTLE with the values as they are in the data file
# intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                  inits = NULL, # Initial parameters = 0
#                                  file = NULL, # Don't save
#                                  # debug = 1, # 1 = estimate, 0 = don't estimate
#                                  random_rec = FALSE, # No random recruitment
#                                  msmMode = 1, # Multispecies mode
#                                  phase = "default")


# Run CEATTLE with differing diet weight proportions --------------------------
# Set different diet weight proportion distributions
wt05 <- c(0.0, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005)
wt10 <- c(0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
wt50 <- c(0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
wt80 <- c(0.0, 0.16, 0.24, 0.32, 0.40, 0.48, 0.56, 0.64, 0.72, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)

# Plot stomach contents curves
# Pull out data from base intrasp run
wts <- hake_intrasp$UobsWtAge %>% 
  group_by(Pred_age) %>%
  summarize(wt_prop = mean(Stomach_proportion_by_weight))

prop <- as.data.frame(cbind(1:15, wt05, wt10, wt50, wt80, wts$wt_prop))
colnames(prop)[c(1, 6)] <- c("age", "empirical data")
prop_all <- melt(prop, id.vars = "age")

stomach_props <- ggplot(prop_all, aes(x=age, y=value, fill=variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  ylab("stomach proportion")
stomach_props

# ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/sensitivity_prop.png", stomach_props,
#        bg = "transparent", width=150, height=80, units="mm", dpi=300)


# Adapt weight proportions to replace those in the excel file & run CEATTLE
run_ceattle <- function(wt, df) {
  new_wt <- c()
  for(i in wt) {
    new_wt <- append(new_wt, c(i, rep(0, 14)))
  }
  df$UobsWtAge$Stomach_proportion_by_weight <- new_wt
  ceattle <- Rceattle::fit_mod(data_list = df,
                               inits = NULL, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # Multi-species mode
                               phase = "default")
  return(ceattle)
}

# Run low-cannibalism models
run_wt05 <- run_ceattle(wt05, hake_intrasp)
run_wt10 <- run_ceattle(wt10, hake_intrasp)
run_wt50 <- run_ceattle(wt50, hake_intrasp)
run_wt80 <- run_ceattle(wt80, hake_intrasp)


# Plot biomass in comparison to no diet & asssessment -------------------------
years <- 1988:2022

# Pull out SSB & overall biomass from CEATTLE runs
ceattle_biomass <- function(run, name) {
  ssb <- (c(run$quantities$biomassSSB) * 2)
  biom <- c(run$quantities$biomass)
  wide <- as.data.frame(cbind(years, ssb, biom))
  colnames(wide) <- c("year", "SSB", "Total Biomass")
  all_biom <- melt(wide, id.vars = "year")
  colnames(all_biom)[2:3] <- c("type", name)
  
  return(all_biom)
}

test_biom <- cbind(ceattle_biomass(run_wt05, "CEATTLE - 0.5% cannibalism"),
                   ceattle_biomass(run_wt10, "CEATTLE - 10% cannibalism"),
                   ceattle_biomass(run_wt50, "CEATTLE - 50% cannibalism"),
                   ceattle_biomass(run_wt80, "CEATTLE - 80% cannibalism"))
test_biom <- test_biom[, c(1:3, 6, 9, 12)]

# Read in intra-species predation run data
intrasp_biom <- read.csv("data/ceattle_intrasp_biomass.csv")
colnames(intrasp_biom)[3] <- "CEATTLE - diet data"

intrasp_R <- read.csv("data/ceattle_intrasp_R.csv")

test_plot_popdy <- function() {
  # Reshape biomass data
  wide <- cbind(test_biom, intrasp_biom[, 3])
  colnames(wide)[ncol(wide)] <- c("CEATTLE - diet data")
  biom <- melt(wide, id.vars = c("year", "type"))
  
  # Put recruitment together
  R_test_all <- cbind(c(run_wt05$quantities$R), c(run_wt10$quantities$R),  
                      c(run_wt50$quantities$R), c(run_wt80$quantities$R))
  R_test_wide <- as.data.frame(cbind(years, R_test_all, intrasp_R))
  colnames(R_test_wide) <- c("year",
                             "CEATTLE - 0.5% cannibalism",
                             "CEATTLE - 10% cannibalism",
                             "CEATTLE - 50% cannibalism",
                             "CEATTLE - 80% cannibalism",
                             "CEATTLE - diet data")
  R_test <- melt(R_test_wide, id.vars = "year")
  
  R_new <- cbind(year = R_test$year,
                 type = rep("Recruitment"),
                 variable = as.character(R_test$variable),
                 value = R_test$value)
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(biom, R_new)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value)

  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = variable, fill = variable)) +
    geom_line(aes(linetype = variable)) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed"), name = "model") +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 1, scales = "free_y")
  
  return(popdy_plot)
}

test_popdy_plot <- test_plot_popdy()
test_popdy_plot

ggsave(filename="plots/CEATTLE/intraspecies predation/Testing/test_intrasp_popdy.png", test_popdy_plot, 
       bg = "transparent", width=170, height=140, units="mm", dpi=300)


# Numbers-at-age for each model run -------------------------------------------
# Read in data from no diet CEATTLE run
intrasp_nbyage <- read.csv("data/ceattle_intrasp_nbyage.csv")
intrasp_nbyage <- cbind(intrasp_nbyage[, -4], rep("CEATTLE - diet data", nrow(intrasp_nbyage)))
colnames(intrasp_nbyage)[4] <- "model"

extract_nbyage <- function(run, name) {
  df <- as.data.frame(as.table(run$quantities$NByage))
  
  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(1988:2022)
  colnames(df) <- c("age", "year", "numbers")
  
  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"
  
  return(df)
}

nbyage_test_all <- rbind(extract_nbyage(run_wt05, "CEATTLE - 0.5% cannibalism"),
                         extract_nbyage(run_wt10, "CEATTLE - 10% cannibalism"),
                         extract_nbyage(run_wt50, "CEATTLE - 50% cannibalism"),
                         extract_nbyage(run_wt80, "CEATTLE - 80% cannibalism"),
                         intrasp_nbyage)

# Set 15 as accumulation age
nbyage_test_all$age[as.numeric(nbyage_test_all$age) > 15] <- 15

# Calculate mean numbers at age & plot
nbyage_test_mean <- nbyage_test_all %>% group_by(age, model) %>%
  summarize(mean_number = mean(numbers))

test_nbyage_plot <- ggplot(nbyage_test_mean, aes(x=age, y=mean_number, fill=model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(labels = c(1:14, "15+")) +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  xlab("age") + ylab("numbers") 
test_nbyage_plot

ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_intrasp_nbyage.png", test_nbyage_plot, 
       bg = "transparent", width=200, height=100, units="mm", dpi=300)


### Compare survey biomass estimate from CEATTLE to true values ---------------
survey <- read.csv("data/assessment/survey_data.csv")
survey <- cbind(survey, model = rep("SS3", length(survey$year)))

extract_srv <- function(run, name){
  df <- data.frame(year = 1995:2019,
                   biomass = run$quantities$srv_bio_hat,
                   log_sd = run$quantities$srv_log_sd_hat,
                   model = rep(name, length(1995:2019)))
  return(df)
}

srv_test <- rbind(extract_srv(run_wt05, "CEATTLE - 0.5% cannibalism"),
                  extract_srv(run_wt10, "CEATTLE - 10% cannibalism"),
                  extract_srv(run_wt50, "CEATTLE - 50% cannibalism"),
                  extract_srv(run_wt80, "CEATTLE - 80% cannibalism"),
                  survey)

test_survey_plot <- ggplot(srv_test, aes(x=year, y=biomass, color=model)) +
  geom_line(alpha = 0.3) +
  geom_point() +
  # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
  scale_color_viridis(discrete = TRUE, direction = -1) +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  xlab("year") + ylab("survey biomass") 
test_survey_plot

ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_intrasp_survey.png", test_survey_plot, 
       bg = "transparent", width=200, height=120, units="mm", dpi=300)


### Look at mortality-at-age timeseries ---------------------------------------
# Custom mortality function that outputs ggplot object ------------------------
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
              p = ggplot2::ggplot(data, aes(y = Age, x = Year, zmin = zlim[1], zmax = zlim[2])) + geom_tile(aes(fill = M))  + scale_y_continuous(expand = c(0, 0), breaks=seq(0,max(ages),round(nages[sp]/5))) + coord_equal() +  scale_x_continuous(expand = c(0, 0))+ theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
              if(!is.null(title)){
                p = p + ggtitle(paste0(title,": ",spnames[j] )) + theme(plot.title = element_text(size = title_cex))
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
m_test <- ggpubr::ggarrange(plot_mortality_custom(Rceattle = run_wt05, type = 0, title = "0.5% cannibalism", maxage = 15),
                            plot_mortality_custom(Rceattle = run_wt10, type = 0, title = "10% cannibalism", maxage = 15),
                            plot_mortality_custom(Rceattle = run_wt50, type = 0, title = "50% cannibalism", maxage = 15),
                            plot_mortality_custom(Rceattle = run_wt80, type = 0, title = "80% cannibalism", maxage = 15),
                            ncol = 2, nrow = 2)

ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_instrasp_M.png", m_test, 
       width=200, height = 100, units = "mm", dpi=300)

