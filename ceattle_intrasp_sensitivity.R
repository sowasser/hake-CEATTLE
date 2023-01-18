# Run CEATTLE for Pacific hake while testing different proportions of
# intraspecies predation

# devtools::install_github("grantdadams/Rceattle@dev")
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230111.xlsx")

# # Run CEATTLE with no predation (single-species)
# nodiet_run <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                 inits = NULL, # Initial parameters = 0
#                                 file = NULL, # Don't save
#                                 msmMode = 0, # Multispecies mode
#                                 phase = "default")


# Run CEATTLE with differing diet weight proportions --------------------------
# Pull out data from base intrasp run
wts <- hake_intrasp$UobsWtAge %>% 
  group_by(Pred_age, Prey_age) %>%
  summarize(wt_prop = mean(Stomach_proportion_by_weight))

wt05 <- rescale_max(wts$wt_prop, to = c(0, 0.005))
wt10 <- rescale_max(wts$wt_prop, to = c(0, 0.1))
wt50 <- rescale_max(wts$wt_prop, to = c(0, 0.5))
wt75 <- rescale_max(wts$wt_prop, to = c(0, 0.75))

prop <- as.data.frame(cbind(wts, wt05 = wt05, wt10 = wt10, wt50 = wt50, wt75 = wt75))
colnames(prop)[3] <- c("observed data")
prop_all <- melt(prop, id.vars = c("Pred_age", "Prey_age"))

# stomach_props <- ggplot(prop_all, aes(x=Prey_age, y=value, fill=variable)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
#   ylab("stomach proportion") + xlab("prey age") +
#   facet_wrap(~Pred_age, ncol = 3)
# stomach_props
# 
# ggsave(filename = "plots/CEATTLE/cannibalism/Testing/sensitivity_prop.png", stomach_props,
#        width=140, height=150, units="mm", dpi=300)


# Adapt weight proportions to replace those in the excel file & run CEATTLE
run_ceattle <- function(wt, df) {
  df$UobsWtAge$Stomach_proportion_by_weight <- wt
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
run_intrasp <- run_ceattle(hake_intrasp$UobsWtAge$Stomach_proportion_by_weight, hake_intrasp)
run_wt05 <- run_ceattle(wt05, hake_intrasp)
run_wt10 <- run_ceattle(wt10, hake_intrasp)
run_wt50 <- run_ceattle(wt50, hake_intrasp)
run_wt75 <- run_ceattle(wt75, hake_intrasp)


# Check fit of CEATTLE model --------------------------------------------------
fit_CEATTLE <- function(run) {
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  jnll_summary <- as.data.frame(run$quantities$jnll_comp)
  jnll_summary$sum <- rowSums(run$quantities$jnll_comp)
  return(list(fit, jnll_summary))
}

sensitivity_fit <- rbind(cbind(model = "wt05", fit_CEATTLE(run_wt05)[[1]]),
                         cbind(model = "wt10", fit_CEATTLE(run_wt10)[[1]]),
                         cbind(model = "wt50", fit_CEATTLE(run_wt50)[[1]]),
                         cbind(model = "wt75", fit_CEATTLE(run_wt75)[[1]]))
sensitivity_summary <- cbind(fit_CEATTLE(run_wt05)[[2]],
                             fit_CEATTLE(run_wt10)[[2]][, 3],
                             fit_CEATTLE(run_wt50)[[2]][, 3],
                             fit_CEATTLE(run_wt75)[[2]][, 3])[, -c(1:2)]
colnames(sensitivity_summary) <- c("wt05", "wt10", "wt50", "wt75")


# Plot biomass & recruitment in comparison to original diet run ---------------
years <- 1988:2019

test_plot_popdy <- function() {
  # Pull out SSB & overall biomass from CEATTLE runs
  ceattle_biomass <- function(run, name) {
    ssb <- (c(run$quantities$biomassSSB) * 2)
    biom <- c(run$quantities$biomass)
    wide <- as.data.frame(cbind(years, ssb, biom))
    colnames(wide) <- c("year", "SSB", "Total Biomass")
    all_biom <- melt(wide, id.vars = "year")
    all_biom2 <- cbind(all_biom,
                       error = c(run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")], 
                                 run$sdrep$sd[which(names(run$sdrep$value) == "biomass")]),
                       model = rep(name))
    return(all_biom2)
  }
  
  test_biom <- rbind(ceattle_biomass(run_intrasp, "observed proportion"),
                     ceattle_biomass(run_wt05, "0.5% cannibalism"),
                     ceattle_biomass(run_wt10, "10% cannibalism"),
                     ceattle_biomass(run_wt50, "50% cannibalism"),
                     ceattle_biomass(run_wt75, "75% cannibalism"))
  
  # Put recruitment together
  ceattle_R <- function(run, name) {
    R <- c(run$quantities$R)
    error <- c(run$sdrep$sd[which(names(run$sdrep$value) == "R")])
    R_all <- as.data.frame(cbind(year = years, 
                                 variable = rep("Recruitment"),
                                 value = R, 
                                 error = error, 
                                 model = rep(name)))
    return(R_all)
  }
  R_test <- rbind(ceattle_R(run_intrasp, "observed proportion"),
                  ceattle_R(run_wt05, "0.5% cannibalism"),
                  ceattle_R(run_wt10, "10% cannibalism"),
                  ceattle_R(run_wt50, "50% cannibalism"),
                  ceattle_R(run_wt75, "75% cannibalism"))
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(test_biom, R_test)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value) / 1000000  # to mt/millions
  all_popdy$error <- as.numeric(all_popdy$error) / 1000000  # to mt/millions
  all_popdy$variable <- factor(all_popdy$variable, labels = c("SSB (mt)", "Total Biomass (mt)", "Recruitment (millions)"))

  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line(aes(linetype = model)) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed"), name = "model") +
    geom_ribbon(aes(ymin=(value-(2*error)), ymax=(value+(2*error))), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")
  
  # Plot ratio of SSB:Biomass to look for skewness in age composition
  ratio <- function(run, name) {
    df <- as.data.frame(cbind(year = years,
                              value = (c(run$quantities$biomassSSB) * 2) / c(run$quantities$biomass),
                              model = rep(name)))
    return(df)
  }
  
  ratio_all <- rbind(ratio(run_intrasp, "observed proportion"),
                     ratio(run_wt05, "0.5% cannibalism"),
                     ratio(run_wt10, "10% cannibalism"),
                     ratio(run_wt50, "50% cannibalism"),
                     ratio(run_wt75, "75% cannibalism"))
  ratio_all$year <- as.numeric(ratio_all$year)
  ratio_all$value <- as.numeric(ratio_all$value)

  ratio_plot <- ggplot(ratio_all, aes(x=year, y=value, color=model)) +
    geom_line(aes(linetype = model)) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed"), name = "model") +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    ylab("SSB/Biomass")
  ratio_plot
  
  return(list(all_popdy, popdy_plot, ratio_plot))
}

test_popdy <- test_plot_popdy()
test_popdy[[2]]

ggsave(filename="plots/CEATTLE/cannibalism/Testing/sensitivity_popdy.png", test_popdy[[2]], 
       width=140, height=150, units="mm", dpi=300)

test_popdy[[3]]
ggsave(filename="plots/CEATTLE/cannibalism/Testing/sensitivity_ratio.png", test_popdy[[3]], 
       width=150, height=80, units="mm", dpi=300)


# Numbers-at-age for each model run -------------------------------------------
# Read in data from no diet CEATTLE run
intrasp_nbyage <- read.csv("data/ceattle_intrasp_nbyage.csv")
intrasp_nbyage <- cbind(intrasp_nbyage[, -4], rep("observed proportion", nrow(intrasp_nbyage)))
colnames(intrasp_nbyage)[4] <- "model"

extract_nbyage <- function(run, name) {
  df <- as.data.frame(as.table(run$quantities$NByage))

  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(years)
  colnames(df) <- c("age", "year", "numbers")

  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"

  return(df)
}

nbyage_test_all <- rbind(extract_nbyage(run_wt05, "0.5% cannibalism"),
                         extract_nbyage(run_wt10, "10% cannibalism"),
                         extract_nbyage(run_wt50, "50% cannibalism"),
                         extract_nbyage(run_wt75, "75% cannibalism"),
                         intrasp_nbyage)

# Set 15 as accumulation age
nbyage_test_all$age[as.numeric(nbyage_test_all$age) > 15] <- 15

# Plot yearly nbyage
nbyage_test_all$age <- as.numeric(nbyage_test_all$age)
# nbyage_test_all$year <- as.numeric(nbyage_test_all$year)

test_nbyage_plot <- ggplot(nbyage_test_all, aes(x=year, y=age)) +
  geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
  scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_y_continuous(breaks = seq(1, 15, 2), labels = c(seq(1, 13, 2), "15+")) +
  scale_x_discrete(breaks = seq(1988, 2019, 3)) +
  xlab(" ") + ylab("Age") +
  theme(legend.position = "none") +
  facet_wrap(~model, ncol = 2, scales = "free_x")
test_nbyage_plot

ggsave(filename = "plots/CEATTLE/cannibalism/Testing/sensitivity_nbyage.png", test_nbyage_plot,
       width=220, height=210, units="mm", dpi=300)


### New plot of popdy and numbers-at-age --------------------------------------
# Calculate annual mean age
mean_nbyage_test <- nbyage_test_all %>%
  group_by(year, model) %>%
  summarize(value = weighted.mean(age, numbers)) %>%
  ungroup()

# Add extra columns, reorder, and combine with popdy dataframe
mean_nbyage_test$variable <- rep("Mean Age")
mean_nbyage_test$error <- rep(0)
mean_nbyage_test <- mean_nbyage_test[, c(1, 4, 3, 5, 2)]

popdy_meanage <- rbind(test_popdy[[1]], mean_nbyage_test) %>%
  filter(variable != "SSB (mt)")
popdy_meanage$year <- as.integer(popdy_meanage$year)
popdy_meanage$model <- factor(popdy_meanage$model)

meanage_popdy_plot <- ggplot(popdy_meanage, aes(x=year, y=value, color = model, fill = model)) +
  geom_line(aes(linetype = model)) +
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed"), name = "model") +
  geom_ribbon(aes(ymin=(value-(2*error)), ymax=(value+(2*error))), alpha = 0.2, color = NA) + 
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  ylab(" ") +
  labs(color = "model") +
  facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
  theme(strip.background = element_blank(), strip.placement = "outside")
meanage_popdy_plot

ggsave(filename="plots/CEATTLE/cannibalism/Testing/sensitivity_meanage_popdy.png", meanage_popdy_plot, 
       width=140, height=150, units="mm", dpi=300)


# ### Compare survey biomass estimate from CEATTLE to true values ---------------
# survey <- read.csv("data/assessment/survey_data.csv")
# survey <- cbind(survey, model = rep("SS3", length(survey$year)))
# 
# extract_srv <- function(run, name){
#   df <- data.frame(year = 1995:2019,
#                    biomass = run$quantities$srv_bio_hat,
#                    log_sd = run$quantities$srv_log_sd_hat,
#                    model = rep(name, length(1995:2019)))
#   return(df)
# }
# 
# srv_test <- rbind(extract_srv(run_wt05, "CEATTLE - 0.5% cannibalism"),
#                   extract_srv(run_wt10, "CEATTLE - 10% cannibalism"),
#                   extract_srv(run_wt50, "CEATTLE - 50% cannibalism"),
#                   extract_srv(run_wt80, "CEATTLE - 80% cannibalism"),
#                   survey)
# 
# test_survey_plot <- ggplot(srv_test, aes(x=year, y=biomass, color=model)) +
#   geom_line(alpha = 0.3) +
#   geom_point() +
#   # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
#   scale_color_viridis(discrete = TRUE, direction = -1) +
#   scale_fill_viridis(discrete = TRUE, direction = -1) +
#   xlab("year") + ylab("survey biomass") 
# test_survey_plot
# 
# ggsave(filename = "plots/CEATTLE/intraspecies predation/Testing/test_intrasp_survey.png", test_survey_plot, 
#        bg = "transparent", width=200, height=120, units="mm", dpi=300)


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
                            plot_mortality_custom(Rceattle = run_wt75, type = 0, title = "75% cannibalism", maxage = 15),
                            ncol = 2, nrow = 2)
# m_test

ggsave(filename = "plots/CEATTLE/cannibalism/Testing/sensitivity_M.png", m_test, 
       width=200, height = 100, units = "mm", dpi=300)
