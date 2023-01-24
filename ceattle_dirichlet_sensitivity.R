# Run CEATTLE for Pacific hake while testing different proportions of
# intraspecies predation, corrected with a Dirichlet multinomial distribution.

# devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230111.xlsx")

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
run_ceattle <- function(df, start, end, proj, init) {
  hake_intrasp$UobsWtAge <- df
  hake_intrasp$styr <- start
  hake_intrasp$endyr <- end
  hake_intrasp$projyr <- proj
  ceattle <- Rceattle::fit_mod(data_list = hake_intrasp,
                               inits = init, # Initial parameters = 0
                               file = NULL, # Don't save
                               # debug = 1, # 1 = estimate, 0 = don't estimate
                               random_rec = FALSE, # No random recruitment
                               msmMode = 1, # Multi-species mode
                               phase = "default")
  return(ceattle)
}

# Run low-cannibalism models
run_all <- run_ceattle(hake_intrasp$UobsWtAge, 1988, 2019, 2019, init = NULL)
run_90s <- run_ceattle(dirichlet_90s, 1988, 1999, 1999, init = NULL)
run_recent <- run_ceattle(dirichlet_recent, 2005, 2019, 2019, init = NULL)

# # Run CEATTLE with no diet
# nodiet_run <- Rceattle::fit_mod(data_list = hake_intrasp,
#                                 inits = NULL, # Initial parameters = 0
#                                 file = NULL, # Don't save
#                                 # debug = 1, # 1 = estimate, 0 = don't estimate
#                                 random_rec = FALSE, # No random recruitment
#                                 msmMode = 0, # Single-species mode
#                                 phase = "default")


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

dirichlet_fit <- rbind(cbind(model = "90s", fit_CEATTLE(run_90s)[[1]]),
                       cbind(model = "recent", fit_CEATTLE(run_recent)[[1]]))
dirichlet_summary <- cbind(fit_CEATTLE(run_90s)[[2]],
                           fit_CEATTLE(run_recent)[[2]][, 3])[, -c(1:2)]
colnames(dirichlet_summary) <- c("90s", "recent")


### Plot population dynamics --------------------------------------------------
timing_plot_popdy <- function() {
  # Pull out SSB & overall biomass from CEATTLE runs
  ceattle_biomass <- function(run, name, years) {
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
  
  test_biom <- rbind(ceattle_biomass(run_all, "all years", 1988:2019),
                     ceattle_biomass(run_90s, "1988-1999", 1988:1999),
                     ceattle_biomass(run_recent, "2005-2019", 2005:2019))
  
  # Put recruitment together
  ceattle_R <- function(run, name, years) {
    R <- c(run$quantities$R)
    error <- c(run$sdrep$sd[which(names(run$sdrep$value) == "R")])
    R_all <- as.data.frame(cbind(year = years, 
                                 variable = rep("Recruitment"),
                                 value = R, 
                                 error = error, 
                                 model = rep(name)))
    return(R_all)
  }
  R_test <- rbind(ceattle_R(run_all, "all years", 1988:2019),
                  ceattle_R(run_90s, "1988-1999", 1988:1999),
                  ceattle_R(run_recent, "2005-2019", 2005:2019))
  
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(test_biom, R_test)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value) / 1000000  # to mt/millions
  all_popdy$error <- as.numeric(all_popdy$error) / 1000000  # to mt/millions
  
  # Find mean difference between the model runs
  mean_SEM <- function(model1, model2, stat, years) {
    df1 <- all_popdy %>% 
      filter(model == model1) %>% filter(variable == stat) %>% filter(year %in% years)
    df2 <- all_popdy %>%
      filter(model == model2) %>% filter(variable == stat) %>% filter(year %in% years)
    mean_out <- mean((df1$value) - (df2$value))
    SEM <- sd((df1$value) - (df2$value)) / sqrt(length(df2$value))
    return(c(paste0(model1, " - ", model2, ", ", stat), mean_out, SEM))
  }
  
  mean_SEM_all <- rbind(mean_SEM("1988-1999", "all years", "SSB", 1988:1999),
                        mean_SEM("1988-1999", "all years", "Total Biomass", 1988:1999),
                        mean_SEM("1988-1999", "all years", "Recruitment", 1988:1999),
                        mean_SEM("2005-2019", "all years", "SSB", 2005:2019),
                        mean_SEM("2005-2019", "all years", "Total Biomass", 2005:2019),
                        mean_SEM("2005-2019", "all years", "Recruitment", 2005:2019))
  
  # Plot population dynamics
  all_popdy$variable <- factor(all_popdy$variable, labels = c("SSB (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)
  
  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line(aes(linetype = model)) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed"), name = "model") +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    ylim(0, NA) +
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")
  
  # Plot ratio of SSB:Biomass to look for skewness in age composition
  ratio <- function(run, name, years) {
    df <- as.data.frame(cbind(year = years,
                              value = (c(run$quantities$biomassSSB) * 2) / c(run$quantities$biomass),
                              model = rep(name)))
    return(df)
  }
  
  ratio_all <- rbind(ratio(run_all, "all years", 1988:2019),
                     ratio(run_90s, "1988-1999", 1988:1999),
                     ratio(run_recent, "2005-2019", 2005:2019))
  ratio_all$year <- as.numeric(ratio_all$year)
  ratio_all$value <- as.numeric(ratio_all$value)

  ratio_plot <- ggplot(ratio_all, aes(x=year, y=value, color=model)) +
    geom_line(aes(linetype = model)) +
    scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed"), name = "model") +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    ylab("SSB/Biomass")

  return(list(mean_SEM_all, popdy_plot, ratio_plot))
}

timing_popdy <- timing_plot_popdy()
mean_SEM <- timing_popdy[[1]]

timing_popdy[[2]]

ggsave(filename="plots/CEATTLE/cannibalism/Testing/dirichlet_popdy.png", timing_popdy[[2]], 
       width=140, height=150, units="mm", dpi=300)

timing_popdy[[3]]
ggsave(filename="plots/CEATTLE/cannibalism/Testing/dirichlet_ratio.png", timing_popdy[[3]], 
       width=150, height=80, units="mm", dpi=300)


# Numbers-at-age for each model run -------------------------------------------
# Read in data from no diet CEATTLE run
extract_nbyage <- function(run, name, years) {
  df <- as.data.frame(as.table(run$quantities$NByage))
  
  df <- df[-seq(0, nrow(df), 2), -c(1:2)]
  levels(df$Var3) <- c(1:20)
  levels(df$Var4) <- c(years)
  colnames(df) <- c("age", "year", "numbers")
  
  df <- cbind(df, rep(name, nrow(df)))
  colnames(df)[4] <- "model"
  # years <- as.character(years)
  # df <- df %>% filter(year %in% years)
  
  return(df)
}

nbyage_test_all <- rbind(extract_nbyage(run_90s, "1988-1999", 1988:1999),
                         extract_nbyage(run_all, "all years", 1988:2019),
                         extract_nbyage(run_recent, "2005-2019", 2005:2019))

# Set 15 as accumulation age
nbyage_test_all$age[as.numeric(nbyage_test_all$age) > 15] <- 15

# Plot yearly nbyage
nbyage_test_all$age <- as.numeric(nbyage_test_all$age)
nbyage_test_all$model <- factor(nbyage_test_all$model, levels = c("1988-1999", "all years", "2005-2019"))

test_nbyage_plot <- ggplot(nbyage_test_all, aes(x=year, y=age)) +
  geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
  scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_y_continuous(breaks = seq(1, 15, 2), labels = c(seq(1, 13, 2), "15+")) +
  scale_x_discrete(breaks = seq(1988, 2019, 3)) +
  xlab(" ") + ylab("Age") +
  theme(legend.position = "none") +
  facet_wrap(~model, ncol = 1)
test_nbyage_plot

ggsave(filename = "plots/CEATTLE/cannibalism/Testing/dirichlet_nbyage.png", test_nbyage_plot,
       width=250, height=150, units="mm", dpi=300)


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
              # scale_fill_viridis(limits=c(0, 2)) +  # set upper limit so plots match
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
m_dirichlet <- gridExtra::grid.arrange(plot_mortality_custom(Rceattle = run_all, type = 0, title = "all years", maxage = 15, zlim = c(0,2.2)),
                                       plot_mortality_custom(Rceattle = run_90s, type = 0, title = "1988-1999", maxage = 15, zlim = c(0,2.2)) +
                                         scale_x_continuous(expand = c(0, 0), breaks = c(1988, 1992, 1996)),
                                       plot_mortality_custom(Rceattle = run_recent, type = 0, title = "2005-2019", maxage = 15, zlim = c(0,2.2)),
                                       ncol = 2, nrow = 2, 
                                       layout_matrix = rbind(c(1,1), c(2,3)))

# # Checking limits of the z-axis, sneakily.
# no_zlim <- gridExtra::grid.arrange(plot_mortality_custom(Rceattle = run_all, type = 0, title = "all years", maxage = 15),
#                                    plot_mortality_custom(Rceattle = run_90s, type = 0, title = "1988-1999", maxage = 15) +
#                                      scale_x_continuous(expand = c(0, 0), breaks = c(1988, 1992, 1996)),
#                                    plot_mortality_custom(Rceattle = run_recent, type = 0, title = "2005-2019", maxage = 15, zlim = c(0, 0.19)),
#                                    ncol = 2, nrow = 2, layout_matrix = rbind(c(1,1), c(2,3)))

ggsave(filename = "plots/CEATTLE/cannibalism/Testing/dirichlet_M.png", 
       m_dirichlet, width=180, height = 180, units = "mm", dpi=300)
