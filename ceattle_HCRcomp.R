#' Comparison of CEATTLE with and without an HCR for the entire projection 
#' period (out to 2100).

library(Rceattle)
library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(ggview)
library(ggsidekick)
# Set ggplot theme
theme_set(theme_sleek())

# ### Run and fit the CEATTLE model ---------------------------------------------
# # Read in CEATTLE data from the excel file
# hake_intrasp <- Rceattle::read_data(file = "data/hake_intrasp_230912.xlsx")
# run_CEATTLE <- function(data, M1, prior, init, msm, estMode) {
#   data$est_M1 <- M1  
#   data$projyr <- 2100
#   # data$endyr <- 2019
#   run <- fit_mod(data_list = data,
#                  inits = init,
#                  file = NULL, # Don't save
#                  msmMode = msm, # Single-species mode - no predation mortality
#                  M1Fun = Rceattle::build_M1(M1_model = M1,
#                                             updateM1 = TRUE,
#                                             M1_use_prior = prior,
#                                             M1_prior_mean = 0.2,
#                                             M1_prior_sd = .1),
#                  # proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
#                  estimateMode = estMode,  # 0 = Fit the hindcast model and projection with HCR specified via HCR; 1 = hindcast only
#                  HCR = Rceattle::build_hcr(HCR = 0),
#                  phase = "default",
#                  initMode = 1,
#                  # random_rec = TRUE,
#                  # random_sel = TRUE,
#                  projection_uncertainty = TRUE,
#                  loopnum = 7) 
#   return(run)
# }
# 
# ss_noHCR <- run_CEATTLE(data = hake_intrasp, 
#                         M1 = 1, 
#                         prior = TRUE, 
#                         init = NULL, 
#                         msm = 0, 
#                         estMode = 0)
# save(ss_model, file = "models/ss_noHCR.Rdata")
# 
# ms_noHCR <- run_CEATTLE(data = hake_intrasp, 
#                         M1 = 1, 
#                         prior = TRUE, 
#                         init = ss_model$model$initial_params, 
#                         msm = 1, 
#                         estMode = 0)
# save(ms_model, file = "models/ms_noHCR.Rdata")

### Load in already-run models ------------------------------------------------
load("models/ss_priorM1.Rdata")
load("models/ms_priorM1.Rdata")
load("models/ss_noHCR.Rdata")
load("models/ms_noHCR.Rdata")

years <- ms_priorM1$model$data_list$styr:ms_priorM1$model$data_list$projyr
ceattle_biomass <- function(run, name, HCR) {
  ssb <- c(run$quantities$biomassSSB * 2)
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
  all_biom2$HCR <- HCR
  return(all_biom2)
}

biomass <- rbind.data.frame(ceattle_biomass(ss_priorM1$model, "Single-species", "40-10"),
                            ceattle_biomass(ms_priorM1$model, "Cannibalism", "40-10"),
                            ceattle_biomass(ss_model$model, "Single-species", "None"),
                            ceattle_biomass(ms_model$model, "Cannibalism", "None"))

ceattle_rec <- function(run, name, HCR) {
  rec <- c(run$quantities$R)
  error <- c(run$sdrep$sd[which(names(run$sdrep$value) == "R")])
  rec_out <- cbind.data.frame(year = years,
                           type = "Recruitment",
                           value = rec,
                           error = error,
                           model = name,
                           HCR = HCR)
  return(rec_out)
}

recruitment <- rbind.data.frame(ceattle_rec(ss_priorM1$model, "Single-species", "40-10"),
                                ceattle_rec(ms_priorM1$model, "Cannibalism", "40-10"),
                                ceattle_rec(ss_model$model, "Single-species", "None"),
                                ceattle_rec(ms_model$model, "Cannibalism", "None"))

all_popdy <- rbind.data.frame(biomass, recruitment)

# Add bounds for error & set 0 as minimum for plotting
all_popdy$min <- all_popdy$value - (2 * all_popdy$error)
all_popdy$min[all_popdy$min < 0] <- 0
all_popdy$max <- all_popdy$value + (2 * all_popdy$error)

all_popdy$type <- factor(all_popdy$type, 
                         labels = c("Spawning Biomass (Mt)", "Total Biomass (Mt)", "Recruitment (millions)"))
all_popdy$model <- factor(all_popdy$model, levels = c("Single-species", "Cannibalism"))

popdy_plot <- ggplot(all_popdy, aes(x=year, y=value, 
                                    color = model, fill = model, linetype = HCR)) +
  geom_line() +
  geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.5) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.5) +
  geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  ylim(0, NA) +
  ylab(" ") + xlab("Year") +
  labs(color = "CEATTLE Model", fill = "CEATTLE Model") +
  facet_wrap(~type, ncol = 1, scales = "free_y", strip.position = "left") +
  theme(strip.background = element_blank(), strip.placement = "outside") 
popdy_plot

ggsave(filename="plots/CEATTLE/cannibalism/HCR_comparison.png", 
       popdy_plot, width=160, height=150, units="mm", dpi=300)