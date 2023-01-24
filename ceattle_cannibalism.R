# Run CEATTLE with intraspecies-predation proportions calculated from diet 
# database going back to 1980.

# devtools::install_github("grantdadams/Rceattle@dev")
# library(Rceattle)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggview)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# Read in CEATTLE data from the excel file
hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230111.xlsx")

# Run and fit the CEATTLE model -----------------------------------------------
run_CEATTLE <- function(data, init, msm) {
  run <- Rceattle::fit_mod(data_list = data,
                           inits = init,
                           file = NULL, # Don't save
                           # debug = 1, # 1 = estimate, 0 = don't estimate
                           random_rec = FALSE, # No random selectivity
                           msmMode = msm, # Single-species mode - no predation mortality
                           phase = "default")
  return(run)
}

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

# Multi-species model using cannibalism data
intrasp_run <- run_CEATTLE(hake_intrasp, init = NULL, msm = 1)
# intrasp_run <- run_CEATTLE(hake_intrasp, init = nodiet_run$estimated_params, msm = 1)
intrasp_fit <- fit_CEATTLE(intrasp_run)[[1]]
intrasp_summary <- fit_CEATTLE(intrasp_run)[[2]]

# No diet (single-species run)
nodiet_init <- run_CEATTLE(hake_intrasp, init = NULL, msm = 0)
nodiet_run <- run_CEATTLE(hake_intrasp, init = nodiet_init$estimated_params, msm = 0)
hake_nodiet <- hake_intrasp
hake_nodiet$est_M1 <- 0
nodiet_fixedM <- run_CEATTLE(hake_nodiet, init = nodiet_init$estimated_params, msm = 0)

nodiet_fit <- rbind(cbind(model = "init = NULL", fit_CEATTLE(nodiet_init)[[1]]),
                    cbind(model = "init = est", fit_CEATTLE(nodiet_run)[[1]]),
                    cbind(model = "fixed M1", fit_CEATTLE(nodiet_fixedM)[[1]]))
nodiet_summary <- fit_CEATTLE(nodiet_run)[[2]]


### Rceattle diagnostic plots -------------------------------------------------
# Rceattle::plot_biomass(intrasp_run, add_ci = TRUE)
# Rceattle::plot_index(intrasp_run)
# Rceattle::plot_catch(intrasp_run)
# Rceattle::plot_selectivity(intrasp_run)
# Rceattle::plot_mortality(intrasp_run)
# Rceattle::plot_indexresidual(intrasp_run)
# Rceattle::plot_logindex(intrasp_run)
# Rceattle::plot_recruitment(intrasp_run, add_ci = TRUE)
# Rceattle::plot_comp(intrasp_run)
# Rceattle::plot_srv_comp(intrasp_run)


### Plot biomass & recruitment in comparison to no diet & assessment ----------
years <- 1988:2019
# Pull out SSB & overall biomass from CEATTLE runs
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

biomass <- ceattle_biomass(intrasp_run, "CEATTLE - cannibalism")
nodiet_biomass <- ceattle_biomass(nodiet_run, "CEATTLE - single-species")
recruitment <- c(intrasp_run$quantities$R)

# Pull out SSB & total biomass from stock synthesis & combine, remove pre-1980
ss3_ssb <- cbind(read.table("data/assessment/ssb.txt")[23:54, 2:3], type = rep("SSB"))
ss3_biomass <- cbind(read.table("data/assessment/biomass.txt")[23:54, 2:3], type = rep("Total Biomass"))
ss3_biom <- as.data.frame(cbind(year = rep(years, 2), 
                                rbind(ss3_ssb, ss3_biomass),
                                model = rep("Assessment")))
colnames(ss3_biom)[2:3] <- c("value", "error")
ss3_biom <- ss3_biom[, c(1, 4, 2, 3, 5)]

# Pull out recruitment
nodiet_R <- c(nodiet_run$quantities$R)
ss3_R <- read.table("data/assessment/recruitment.txt")[23:57,]

# Find mean difference between the model runs
mean_SEM <- function(df1, df2, title) {
  mean_out <- mean((df1 / 1000000) - (df2 / 1000000))
  SEM <- sd((df1 / 1000000) - (df2 / 1000000)) / sqrt(length(range))
  return(c(title, mean_out, SEM))
}

mean_SEM_all <- rbind(mean_SEM(biomass[1:32, 4], nodiet_biomass[1:32, 4], "cannibalism - no diet, SSB"),
                      mean_SEM(biomass[33:64, 4], nodiet_biomass[33:64, 4], "cannibalism - no diet, Total"),
                      mean_SEM(recruitment, nodiet_R, "cannibalism - no diet, R"),
                      mean_SEM(nodiet_biomass[1:32, 4], ss3_biom[1:32, 4], "no diet - SS3, SSB"),
                      mean_SEM(nodiet_biomass[33:64, 4], ss3_biom[33:64, 4], "no diet - SS3, Total"),
                      mean_SEM(nodiet_R, ss3_R[1:32, 2], "no diet - SS3, R"))


# Read in other model runs for comparison and plot
plot_popdy <- function() {
  # Put biomass together
  nodiet_biom <- ceattle_biomass(nodiet_run, "CEATTLE - single-species")
  
  # Combine all biomass sources together
  biom_all <- rbind(biomass, nodiet_biom, ss3_biom)
  biom_all$value <- biom_all$value / 1000000  # to mt
  biom_all$error <- biom_all$error / 1000000  # to mt

  # Put recruitment together
  R_wide <- data.frame(year = years, recruitment, nodiet_R)
  colnames(R_wide)[2:3] <- c("CEATTLE - cannibalism", "CEATTLE - single-species")
  R <- melt(R_wide, id.vars = "year")
  # Offset the stock synthesis data by one year (min age in CEATTLE is 1; in SS3 is 0)
  ss3_1 <- as.data.frame(cbind(year = 1989:2019, 
                               variable = rep("Assessment"), 
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
                 value = R_all$value / 1000000,  # to millions
                 error = R_all$error / 1000000,  # to millions
                 model = as.character(R_all$variable))
  
  # Combine biomass & recruitment and plot
  all_popdy <- rbind(biom_all, R_new)
  all_popdy$year <- as.numeric(all_popdy$year)
  all_popdy$value <- as.numeric(all_popdy$value)
  all_popdy$error <- as.numeric(all_popdy$error)
  all_popdy$model <- factor(all_popdy$model, levels = c("Assessment", "CEATTLE - single-species", "CEATTLE - cannibalism"))
  all_popdy$type <- factor(all_popdy$type, labels = c("SSB (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))
  
  # Add bounds for error & set 0 as minimum for plotting
  all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
  all_popdy$min[all_popdy$min < 0] <- 0
  all_popdy$max <- all_popdy$value + (2 * all_popdy$error)

  popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, color = model, fill = model)) +
    geom_line() +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    ylim(0, NA) +
    ylab(" ") +
    labs(color = "model") +
    facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")
  
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
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    ylab("SSB/Biomass")
  
  # Plot the difference between model runs
  ceattle_intrasp <- all_popdy %>% filter(model == "CEATTLE - cannibalism")
  ceattle_nodiet <- all_popdy %>% filter(model == "CEATTLE - single-species")
  assessment <- all_popdy %>% filter(model == "Assessment")

  diff_intrasp <- rbind(cbind.data.frame(year = years,
                                         type = ceattle_intrasp$type,
                                         difference = ceattle_intrasp$value - ceattle_nodiet$value,
                                         models = "cannibalism - single-species"),
                        cbind.data.frame(year = years,
                                         type = ceattle_nodiet$type,
                                         difference = ceattle_nodiet$value - assessment$value,
                                         models = "single-species - assessment"))
  
  diff_plot <- ggplot(diff_intrasp, aes(x=year, y=difference, color=models)) +
    geom_line() +
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside")
  
  return(list(popdy_plot, ratio_plot, diff_plot))

    
}

popdy <- plot_popdy()
popdy[[1]]
popdy[[2]]
popdy[[3]]


### Numbers-at-age for each model run -----------------------------------------
# Extract numbers at age for cannibalism model run
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

nbyage <- extract_byage(intrasp_run$quantities$NByage, "CEATTLE - cannibalism", "numbers")

plot_nbyage <- function() {
  # # Read in data from no diet CEATTLE run
  # nbyage_nodiet <- extract_nbyage(nodiet_run, "CEATTLE - single species")

  # Read in data from SS3 & average beginning & middle of the year
  nbyage_ss3_all <- read.csv("data/assessment/nbyage.csv")[c(9:40, 52:83), ]
  colnames(nbyage_ss3_all) <- c("year", "timing", c(0:20))

  nbyage_ss3_wide <- nbyage_ss3_all %>%
    group_by(year) %>%
    summarize_at(vars("0":"20"), mean)

  nbyage_ss3 <- melt(nbyage_ss3_wide[, -2], id.vars = "year")
  nbyage_ss3 <- cbind(nbyage_ss3, rep("Stock Assessment", length(nbyage_ss3$year)))
  colnames(nbyage_ss3)[2:4] <- c("age", "numbers", "model")

  # Combine with nbyage from intrasp run
  nbyage_all <- rbind(nbyage, nbyage_ss3)
  
  # Find mean numbers-at-age for each year
  nbyage_all$age <- as.numeric(nbyage_all$age)
  nbyage_mean <- nbyage_all %>%
    group_by(year, model) %>%
    summarize(mean = weighted.mean(age, numbers))
  nbyage_mean$year <- as.numeric(nbyage_mean$year)
  
  mean_nbyage_plot <- ggplot(nbyage_mean, aes(x = year, y = mean)) +
    geom_line(aes(color = model)) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) 
  mean_nbyage_plot

  # Set 15 as accumulation age
  nbyage_all$age[as.numeric(nbyage_all$age) > 15] <- 15
  
  # Reduce down to millions for plotting
  nbyage_all$numbers <- nbyage_all$numbers / 1000000

  # Plot yearly nbyage
  nbyage_plot <- ggplot(nbyage_all, aes(x=year, y=age)) +
    geom_point(aes(size = numbers, color = numbers, fill = numbers)) +
    scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
    scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
    scale_y_continuous(breaks = seq(1, 15, 2), labels = c(seq(1, 13, 2), "15+")) +
    scale_x_discrete(breaks = seq(1988, 2019, 3)) +
    xlab(" ") + ylab("Age") + labs(fill="millions (n)", size="millions (n)", color="millions (n)") +
    facet_wrap(~model, ncol=1)
}

nbyage_plot <- plot_nbyage()
nbyage_plot


### Biomass-at-age for each model run -----------------------------------------
biombyage <- extract_byage(intrasp_run$quantities$biomassByage, "CEATTLE - cannibalism", "biomass")

# Set 15 as accumulation age
biombyage$age[as.numeric(biombyage$age) > 15] <- 15
biombyage$age <- as.integer(biombyage$age)
  
# Plot yearly biomass by age
biombyage_plot <- ggplot(biombyage, aes(x=year, y=age)) +
  geom_point(aes(size = biomass, color = biomass, fill = biomass)) +
  scale_fill_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_color_viridis(direction = -1, begin = 0.1, end = 0.9) +
  scale_y_continuous(breaks = seq(1, 15, 2), labels = c(seq(1, 13, 2), "15+")) +
  scale_x_discrete(breaks = seq(1988, 2019, 3)) +
  xlab(" ") + ylab("Age") 
biombyage_plot


### Plot realized consumption -------------------------------------------------
b_consumed <- extract_byage(intrasp_run$quantities$B_eaten_as_prey, "CEATTLE - cannibalism", "biomass")[, -4]
b_consumed$age <- as.integer(b_consumed$age)

# Filter to only ages below 6 (the ages of hake consumed)
b_consumed <- b_consumed %>% filter(age < 6)
b_consumed$age <- as.factor(b_consumed$age)

# Plot yearly biomass consumed by age
b_consumed_plot <- ggplot(b_consumed, aes(x=year, y=(biomass / 1000000), fill = age)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(breaks = seq(1988, 2019, 3)) +
  ylab("Biomass consumed (mt)")
b_consumed_plot

# Plot ratio of biomass consumed to approximation of predator biomass (SSB)
# Add ages together
yearly_consumed <- b_consumed %>%
  group_by(year) %>%
  summarize(total_biomass = sum(biomass)) %>%
  ungroup()

ssb <- biomass %>% filter(type == "SSB")

yearly_consumed$total_biomass <- (yearly_consumed$total_biomass / (ssb$value))
yearly_consumed$year <- as.numeric(as.character(yearly_consumed$year))
yearly_b_plot <- ggplot(yearly_consumed, aes(x = year, y = total_biomass)) +
  geom_line() +
  ylab("biomass of prey / SSB")
yearly_b_plot

# plot_b_eaten(intrasp_run)

### Compare survey biomass estimate from CEATTLE to true values ---------------
survey_biom <- function(run, name) {
  srv <- data.frame(year = 1995:2019,
                    biomass = run$quantities$srv_bio_hat,
                    log_sd = run$quantities$srv_log_sd_hat,
                    model = rep(name, length(1995:2019)))
  return(srv)
}

plot_survey <- function() {
  intrasp_srv <- survey_biom(intrasp_run, "CEATTLE - cannibalism")

  nodiet_srv <- survey_biom(nodiet_run, "CEATTLE - no diet")

  survey <- read.csv("data/assessment/survey_data.csv")
  survey <- cbind(survey, model = rep("Stock Synthesis", length(survey$year)))

  survey_all <- rbind(intrasp_srv, nodiet_srv, survey)

  survey_plot <- ggplot(survey_all, aes(x=year, y=biomass, color=model)) +
    geom_line(alpha = 0.3) +
    geom_point() +
    # geom_ribbon(aes(ymin=(biomass-log_sd), ymax=(biomass+log_sd), fill=model)) +  # Including log sd, but values are really small!
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
    xlab("year") + ylab("survey biomass")

  return(survey_plot)
}

survey_plot <- plot_survey()
survey_plot


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
            return(list(p, data))
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

M <- plot_mortality_custom(Rceattle = intrasp_run, type = 0, title = NULL, maxage = 15, zlim = c(0,1.5)) 
M[[1]]  # mortality plot

# Examine mortality data
intrasp_run$quantities$M1  # 0.317
nodiet_run$quantities$M1  # 0.257

M_data <- M[[2]]  # natural mortality data from the model

# Just mortality for age 1
M_age1 <- M_data %>% filter(Age == 1)  
max(M_age1$M)
min(M_age1$M)

# Mean natural mortality across the time-series for predated ages
M_mean <- M_data %>%
  filter(Age < 6) %>%
  group_by(Year) %>%
  summarize(mean_M = mean(M))
max(M_mean$mean_M)
min(M_mean$mean_M)


### Save data & plots (when not experimenting) --------------------------------
# # Data
# write.csv(biomass, "data/ceattle_intrasp_biomass.csv", row.names = FALSE)
# write.csv(nodiet_biomass, "data/ceattle_nodiet_biomass.csv", row.names = FALSE)
# write.csv(recruitment, "data/ceattle_intrasp_R.csv", row.names = FALSE)
# write.csv(nbyage, "data/ceattle_intrasp_nbyage.csv", row.names = FALSE)
# 
# # Plots
# ggsave(filename="plots/CEATTLE/cannibalism/popdy.png", popdy[[1]], width=140, height=150, units="mm", dpi=300)
# ggsave(filename="plots/CEATTLE/cannibalism/biomass_ratio.png", popdy[[2]], bg = "transparent", width=150, height=80, units="mm", dpi=300)
# ggsave(filename = "plots/CEATTLE/cannibalism/nbyage.png", nbyage_plot, bg = "white", width=160, height=120, units="mm", dpi=300)
# ggsave(filename = "plots/CEATTLE/cannibalism/biomass_byage.png", biombyage_plot, bg = "white", width=160, height=80, units="mm", dpi=300)
# ggsave(filename = "plots/CEATTLE/cannibalism/realized_consumption.png", yearly_b_plot, bg = "white", width=140, height=80, units="mm", dpi=300)
# ggsave(filename = "plots/CEATTLE/cannibalism/survey_biomass.png", survey_plot, bg = "white", width=200, height=120, units="mm", dpi=300)
# ggsave(filename = "plots/CEATTLE/cannibalism/M.png", M[[1]], width = 160, height = 70, units = "mm", dpi=300)
