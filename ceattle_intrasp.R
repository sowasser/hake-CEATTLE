# Run CEATTLE with intraspecies-predation proportions calculated from diet 
# database going back to 1980.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
# Set transparent ggplot theme
source("~/Desktop/Local/ggsidekick/R/theme_sleek_transparent.R")
theme_set(theme_sleek_transparent())
halloween <- c("darkorchid3", "darkorange", "chartreuse3", "deepskyblue3")

# Read in CEATTLE data from the excel file
hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_221026.xlsx")

intrasp_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                                 inits = NULL, # Initial parameters = 0
                                 file = NULL, # Don't save
                                 # debug = 1, # 1 = estimate, 0 = don't estimate
                                 random_rec = FALSE, # No random recruitment
                                 msmMode = 1, # Multi-species mode
                                 phase = "default")

hake_nodiet <- hake_intrasp
hake_nodiet$est_M1 <- 0  # Use base M1
nodiet_run <- Rceattle::fit_mod(data_list = hake_intrasp,
                                inits = NULL, # Initial parameters = 0
                                file = NULL, # Don't save
                                # debug = 1, # 1 = estimate, 0 = don't estimate
                                random_rec = FALSE, # No random recruitment
                                msmMode = 0, # Single-species mode - no predation mortality
                                phase = "default")
# plot_biomass(nodiet_run, add_ci = TRUE)

# # Rceattle diagnostics ------------------------------------------------------
# plot_biomass(intrasp_run, add_ci = TRUE)
# plot_index(intrasp_run)
# plot_catch(intrasp_run)
# plot_selectivity(intrasp_run)
# plot_mortality(intrasp_run)
# plot_indexresidual(intrasp_run)
# plot_logindex(intrasp_run)
# plot_recruitment(intrasp_run, add_ci = TRUE)
# plot_comp(intrasp_run)
# plot_srv_comp(intrasp_run)


### Plot biomass & recruitment in comparison to no diet & assessment ----------
years <- 1988:2019
# Pull out SSB & overall biomass from CEATTLE runs
ceattle_biomass <- function(run, name) {
  # run <- intrasp_run
  # name <- "CEATTLE - cannibalism"
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

biomass <- ceattle_biomass(intrasp_run, "CEATTLE - cannibalism")
write.csv(biomass, "data/ceattle_intrasp_biomass.csv", row.names = FALSE)

recruitment <- c(intrasp_run$quantities$R)
write.csv(recruitment, "data/ceattle_intrasp_R.csv", row.names = FALSE)


# Read in other model runs for comparison and plot
plot_popdy <- function() {
  # Put biomass together
  nodiet_biom <- ceattle_biomass(nodiet_run, "CEATTLE - single-species")
  # Pull out SSB & total biomass from stock synthesis & combine, remove pre-1980
  ss3_ssb <- cbind(read.table("data/assessment/ssb.txt")[23:54, 2:3], type = rep("SSB"))
  ss3_biomass <- cbind(read.table("data/assessment/biomass.txt")[23:54, 2:3], type = rep("Total Biomass"))
  ss3_biom <- as.data.frame(cbind(year = rep(years, 2), 
                                  rbind(ss3_ssb, ss3_biomass),
                                  model = rep("Assessment", length(ss3_ssb$V2) * 2)))
  colnames(ss3_biom)[2:3] <- c("value", "error")
  ss3_biom <- ss3_biom[, c(1, 4, 2, 3, 5)]
  biom_all <- rbind(biomass, nodiet_biom, ss3_biom)
  
  # Put recruitment together
  nodiet_R <- c(nodiet_run$quantities$R)
  ss3_R <- read.table("data/assessment/recruitment.txt")[23:57,]
  R_wide <- data.frame(year = years, recruitment, nodiet_R)
  colnames(R_wide)[2:3] <- c("CEATTLE - cannibalism", "CEATTLE - single-species")
  R <- melt(R_wide, id.vars = "year")
  # Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS3 is 0)
  ss3_1 <- as.data.frame(cbind(year = 1989:2019, 
                               variable = rep("Assessment", (length(1989:2019))), 
                               value = ss3_R[1:length(1989:2019), 2],
                               error = ss3_R[1:length(1989:2019), 3]))
  R_all <- rbind(cbind(R, error = c(intrasp_run$sdrep$sd[which(names(intrasp_run$sdrep$value) == "R")], 
                                    nodiet_run$sdrep$sd[which(names(nodiet_run$sdrep$value) == "R")])), 
                 ss3_1)
  R_all$value <- as.numeric(R_all$value)
  R_all$year <- as.numeric(R_all$year)
  R_all$error <- as.numeric(R_all$error)
  # Reshape to match biomass
  R_new <- cbind(year = R_all$year, 
                 type = rep("Recruitment"), 
                 value = R_all$value,
                 error = R_all$error,
                 model = as.character(R_all$variable))
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(biom_all, R_new)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value)
  all_popdy$error <- as.numeric(all_popdy$error)
  all_popdy$model <- factor(all_popdy$model, levels = c("Assessment", "CEATTLE - single-species", "CEATTLE - cannibalism"))
  
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line() +
    geom_ribbon(aes(ymin=(value-(2*error)), ymax=(value+(2*error))), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 2, scales = "free_y")
  
  # Plot ratio of SSB:Biomass to look for skewness in age composition
  ratio <- as.data.frame(cbind(year = years,
                               assessment = ss3_ssb[, 1]/ss3_biomass[, 1],
                               cannibalism = biomass[1:length(years), 3] / biomass[(length(years)+1):(length(years)*2), 3],
                               single_species = nodiet_biom[1:length(years), 3] / nodiet_biom[(length(years)+1):(length(years)*2), 3]))
                 
  colnames(ratio) <- c("year", "Assessment", "CEATTLE - cannibalism", "CEATTLE - single-species")
  ratio2 <- melt(ratio, id.vars = "year", variable.name = "model")
  ratio2$model <- factor(ratio2$model, levels = c("Assessment", "CEATTLE - single-species", "CEATTLE - cannibalism"))
  
  ratio_plot <- ggplot(ratio2, aes(x=year, y=value, color=model)) +
    geom_line() +
    scale_color_manual(values = halloween) +
    ylab("SSB/Biomass")
  
  return(list(popdy_plot, ratio_plot))
}

popdy <- plot_popdy()
popdy[[1]]

ggsave(filename="plots/CEATTLE/intraspecies predation/intrasp_popdy.png", popdy[[1]], 
       bg = "transparent", width=280, height=140, units="mm", dpi=300)

popdy[[2]]
ggsave(filename="plots/CEATTLE/intraspecies predation/intrasp_ratio.png", popdy[[2]], 
       bg = "transparent", width=150, height=80, units="mm", dpi=300)


# ### Numbers-at-age for each model run -----------------------------------------
# # Extract numbers at age for intraspecies predation model run
# extract_nbyage <- function(run, name) {
#   df <- as.data.frame(as.table(run$quantities$NByage))
#   
#   df <- df[-seq(0, nrow(df), 2), -c(1:2)]
#   levels(df$Var3) <- c(1:20)
#   levels(df$Var4) <- c(years)
#   colnames(df) <- c("age", "year", "numbers")
#   
#   df <- cbind(df, rep(name, nrow(df)))
#   colnames(df)[4] <- "model"
#   
#   return(df)
# }
# 
# nbyage <- extract_nbyage(intrasp_run, "CEATTLE - cannibalism")
# write.csv(nbyage, "data/ceattle_intrasp_nbyage.csv", row.names = FALSE)
# 
# plot_nbyage <- function(output) {
#   # Read in data from no diet CEATTLE run
#   nbyage_nodiet <- extract_nbyage(nodiet_run, "CEATTLE - single species")
# 
#   # Read in data from SS3 & average beginning & middle of the year
#   nbyage_ss3_all <- read.csv("data/assessment/nbyage.csv")[c(9:40, 52:83), ]
#   colnames(nbyage_ss3_all) <- c("year", "timing", c(0:20))
#   
#   nbyage_ss3_wide <- nbyage_ss3_all %>%
#     group_by(year) %>%
#     summarize_at(vars("0":"20"), mean)
#   
#   nbyage_ss3 <- melt(nbyage_ss3_wide[, -2], id.vars = "year")
#   nbyage_ss3 <- cbind(nbyage_ss3, rep("Stock Synthesis", length(nbyage_ss3$year)))
#   colnames(nbyage_ss3)[2:4] <- c("age", "numbers", "model")
#   
#   # Combine with nbyage from intrasp run
#   nbyage_all <- rbind(nbyage, nbyage_nodiet, nbyage_ss3)
#   
#   # Set 15 as accumulation age
#   nbyage_all$age[as.numeric(nbyage_all$age) > 15] <- 15
#   
#   # Calculate mean numbers at age & plot
#   nbyage_mean <- nbyage_all %>% 
#     group_by(age, model) %>%
#     summarize(mean = mean(numbers), sd = sd(numbers)) 
#   
#   # Plot mean nbyage across years
#   nbyage_plot_mean <- ggplot(nbyage_mean, aes(x=age, y=mean, fill=model, color=model)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2, position=position_dodge(.9)) +  # only upper error bars
#     scale_x_discrete(labels = c(1:14, "15+")) +
#     scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
#     scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
#     xlab("age") + ylab("numbers") 
#   
#   # Plot yearly nbyage
#   nbyage_plot_yearly <- ggplot(nbyage_all, aes(x=age, y=numbers, fill=model)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     scale_x_discrete(labels = c(1:14, "15+")) +
#     scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
#     xlab("age") + ylab("numbers") +
#     facet_wrap(~ year)
#   
#   # Conditionally return plots
#   if(output == "mean") {return(nbyage_plot_mean)}
#   if(output == "yearly") {return(nbyage_plot_yearly)}
# }
# 
# nbyage_plot_mean <- plot_nbyage(output = "mean")
# nbyage_plot_mean
# 
# nbyage_plot_yearly <- plot_nbyage(output = "yearly")
# nbyage_plot_yearly
# 
# ggsave(filename = "plots/CEATTLE/intraspecies predation/nbyage_intrasp.png", nbyage_plot_mean, 
#        bg = "white", width=170, height=90, units="mm", dpi=300)
# 
# 
# ### Compare survey biomass estimate from CEATTLE to true values ---------------
# survey_biom <- function(run, name) {
#   srv <- data.frame(year = 1995:2019,
#                     biomass = run$quantities$srv_bio_hat,
#                     log_sd = run$quantities$srv_log_sd_hat,
#                     model = rep(name, length(1995:2019)))
#   return(srv)
# }
# 
# plot_survey <- function() {
#   intrasp_srv <- survey_biom(intrasp_run, "CEATTLE - cannibalism")
# 
#   nodiet_srv <- survey_biom(nodiet_run, "CEATTLE - no diet")
#   
#   survey <- read.csv("data/assessment/survey_data.csv")
#   survey <- cbind(survey, model = rep("Stock Synthesis", length(survey$year)))
#   
#   survey_all <- rbind(intrasp_srv, nodiet_srv, survey)
#   
#   survey_plot <- ggplot(survey_all, aes(x=year, y=biomass, color=model)) +
#     geom_line(alpha = 0.3) +
#     geom_point() +
#     # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
#     scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
#     scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
#     xlab("year") + ylab("survey biomass") 
#   
#   return(survey_plot)
# }
# 
# survey_plot <- plot_survey()
# survey_plot
# 
# ggsave(filename = "plots/CEATTLE/intraspecies predation/survey_biomass.png", survey_plot, 
#        bg = "white", width=200, height=120, units="mm", dpi=300)


### Compare predation mortality (M2) ------------------------------------------
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

m_plot <- plot_mortality_custom(Rceattle = intrasp_run, type = 0, title = NULL, maxage = 15) 
m_plot

ggsave(filename = "plots/CEATTLE/intraspecies predation/M_intrasp.png", m_plot, 
       width = 160, height = 70, units = "mm", dpi=300)

