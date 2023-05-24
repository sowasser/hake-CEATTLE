# Profile over M1 

# Re-install Rceattle if needed
# remove.packages("Rceattle")
# remove.packages("00LOCK-Rceattle")
# devtools::install_github("grantdadams/Rceattle", ref = "dev")

library(Rceattle)
library(ggplot2)
library(viridis)
library(dplyr)

data <- read_data(file = "data/hake_intrasp_230427.xlsx")  # Read in data

# Function updating M1 for each run of the model
get_profile <- function(M1_change, init, msm) {
  data$M1_base[, 3:17] <- M1_change
  run <- fit_mod(
    data_list = data,
    inits = init,
    file = NULL, # Don't save
    msmMode = 0, 
    M1Fun = Rceattle::build_M1(M1_model = 0,
                               updateM1 = TRUE,
                               M1_use_prior = FALSE),
    estimateMode = 0,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
    HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                              FsprLimit = 0.4, # F40%
                              Ptarget = 0.4, # Target is 40% B0
                              Plimit = 0.1, # No fishing when SB<SB10
                              Pstar = 0.45,
                              Sigma = 0.5),
    phase = "default",
    initMode = 1,
    projection_uncertainty = TRUE,
    getsd = FALSE
  )
  
  # Save the resulting run to losing progress to R bombs!
  if (msm == 0) {save(run, file = paste0("models/profile/ss/run", as.character(M1_change), ".Rdata"))}
  if (msm == 1) {save(run, file = paste0("models/profile/ms/run", as.character(M1_change), ".Rdata"))}
  
  message(run$quantities$jnll)  # print JNLL to console
  return(run)
}

### Run profile over M1 -------------------------------------------------------
# Load model with estimated M1 & check starting value
load("models/ss_estM1.Rdata")
startM_ss <- round(ss_estM1$model$quantities$M1[1, 1, 1], digits = 2)

# SS down
run0 <- get_profile(0.26, ss_estM1$model$estimated_params, 0)
run1 <- get_profile(0.25, run0$estimated_params, 0)
run2 <- get_profile(0.24, run1$estimated_params, 0)
run3 <- get_profile(0.23, run2$estimated_params, 0)
run4 <- get_profile(0.22, run3$estimated_params, 0)
run5 <- get_profile(0.21, run4$estimated_params, 0)
run6 <- get_profile(0.20, run5$estimated_params, 0)
run7 <- get_profile(0.19, run6$estimated_params, 0)
run8 <- get_profile(0.18, run7$estimated_params, 0)
run9 <- get_profile(0.17, run8$estimated_params, 0)
run10 <- get_profile(0.16, run9$estimated_params, 0)
run11 <- get_profile(0.15, run10$estimated_params, 0)

# SS up
run12 <- get_profile(0.27, run$estimated_params, 0)
run13 <- get_profile(0.28, run12$estimated_params, 0)
run14 <- get_profile(0.29, run13$estimated_params, 0)
run15 <- get_profile(0.30, run14$estimated_params, 0)
run16 <- get_profile(0.31, run15$estimated_params, 0)
run17 <- get_profile(0.32, run16$estimated_params, 0)
run18 <- get_profile(0.33, run17$estimated_params, 0)
run19 <- get_profile(0.34, run18$estimated_params, 0)
run20 <- get_profile(0.35, run19$estimated_params, 0)

rm(list = ls())  # clear environment to re-set runs

# Load model with estimated M1 & check starting value
load("models/ms_estM1.Rdata")
startM_ms <- round(ms_estM1$model$quantities$M1[1, 1, 1], digits = 2)

# MS down
run0 <- get_profile(0.32, ms_estM1$model$estimated_params, 1)
run1 <- get_profile(0.31, run0$estimated_params, 1)
run2 <- get_profile(0.30, run1$estimated_params, 1)
run3 <- get_profile(0.29, run2$estimated_params, 1)
run4 <- get_profile(0.28, run3$estimated_params, 1)
run5 <- get_profile(0.27, run4$estimated_params, 1)
run6 <- get_profile(0.26, run5$estimated_params, 1)
run7 <- get_profile(0.25, run6$estimated_params, 1)
run8 <- get_profile(0.24, run7$estimated_params, 1)
run9 <- get_profile(0.23, run8$estimated_params, 1)
run10 <- get_profile(0.22, run9$estimated_params, 1)
run11 <- get_profile(0.21, run10$estimated_params, 1)
run12 <- get_profile(0.20, run11$estimated_params, 1)
run13 <- get_profile(0.19, run12$estimated_params, 1)
run14 <- get_profile(0.18, run13$estimated_params, 1)
run15 <- get_profile(0.17, run14$estimated_params, 1)
run16 <- get_profile(0.16, run15$estimated_params, 1)
run17 <- get_profile(0.15, run16$estimated_params, 1)

# MS up
run18 <- get_profile(0.33, run0$estimated_params, 1)
run19 <- get_profile(0.34, run18$estimated_params, 1)
run20 <- get_profile(0.35, run19$estimated_params, 1)


### Get JNLL for each run and plot --------------------------------------------
# Clear environment and load models back in
rm(list = ls())  
load("models/ss_estM1.Rdata")
load("models/ms_estM1.Rdata")

# Single species
runs_ss <- list.files(path = "models/profile/ss")  # List of all model runs

# Get sum of joint negative log likelihood for each run
jnll_all_ss <- c() 
m1_all_ss <- c()
for(i in 1:length(runs_ss)) {
  load(paste0("models/profile/ss/", runs_ss[i]))
  jnll_all_ss[i] <- run$quantities$jnll
  m1_all_ss[i] <- round(run$quantities$M1[1, 1, 1], digits = 2)
}

# Combine with M1 input for each run
profile_ss <- cbind.data.frame(M1 = m1_all_ss,
                               JNLL = jnll_all_ss - ss_estM1$model$quantities$jnll,
                               model = "single-species")

# Cannibalism
runs_ms <- list.files(path = "models/profile/ms")  # List of all model runs

# Get sum of joint negative log likelihood for each run
jnll_all_ms <- c() 
m1_all_ms <- c()
for(i in 1:length(runs_ms)) {
  load(paste0("models/profile/ms/", runs_ms[i]))
  jnll_all_ms[i] <- run$quantities$jnll
  m1_all_ms[i] <- round(run$quantities$M1[1, 1, 1], digits = 2)
}

# Combine with M1 input for each run
profile_ms <- cbind.data.frame(M1 = seq(0.15, 0.35, 0.01),
                               JNLL = jnll_all_ms - ms_estM1$model$quantities$jnll,
                               model = "cannibalism")

est_M1 <- cbind.data.frame(value = c(0.257, 0.318),
                           model = c("single-species", "cannibalism"))
est_M1$model <- factor(est_M1$model, levels = c("single-species", "cannibalism"))

profile_plot <- rbind(profile_ss, profile_ms) %>%
  mutate(model = factor(model, levels = c("single-species", "cannibalism"))) %>%
  ggplot(., aes(x = M1, y = JNLL)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0.21, linetype = "dotted") +
  geom_vline(data = est_M1, mapping = aes(xintercept = value, color = model), linetype = "dashed") +
  ylab("JNLL") + labs(color = "M1 estimate") +
  # scale_x_continuous(breaks = seq(0.15, 0.35, 0.01)) +  # all scale markers for investigating
  ggsidekick::theme_sleek() +
  facet_wrap(~model)
profile_plot

ggsave(filename="plots/CEATTLE/cannibalism/Testing/M1_profile.png", 
       profile_plot, 
       width=180, height=80, units="mm", dpi=300)

### Plot JNLL components ------------------------------------------------------
comp_out <- function(run) {
  comp <- data.frame(run$quantities$jnll_comp)
  comp$component <- rownames(comp)
  rownames(comp) <- NULL
  comp[nrow(comp) + 1, ] <- c(sum(comp[, 1]), sum(comp[, 2]), "Total NLL")  # add total NLL
  # Separate comps for fishery & survey (in different columns originally)
  comp[nrow(comp) + 1, ] <- c(comp[3, 1], 0, "Fishery age composition")
  comp[nrow(comp) + 1, ] <- c(comp[3, 2], 0, "Survey age composition")
  comp$NLL <- as.numeric(comp$Sp.Srv.Fsh_1) + as.numeric(comp$Sp.Srv.Fsh_2)  # combine species together
  comp <- comp[, c("component", "NLL")]
  comp <- comp %>% filter(NLL != 0)  # remove components w/ no likelihood
  comp$NLL <- as.numeric(comp$NLL)
  # Select components for easier plotting.
  comp <- comp %>% filter(component %in% c("Fishery age composition",
                                           "Survey age composition",
                                           "Recruitment deviates",
                                           "Selectivity deviates",
                                           "Survey biomass",
                                           "Total NLL"))
  return(comp)
}

est_comp_ss <- comp_out(ss_estM1$model)
comp_all_ss <- data.frame()
for(i in 1:length(runs_ss)) {
  load(paste0("models/profile/ss/", runs_ss[i]))
  comp <- comp_out(run)
  comp$NLL <- comp$NLL - est_comp_ss$NLL
  comp$M1 <- round(run$quantities$M1[1, 1, 1], digits = 2)
  comp_all_ss <- rbind(comp_all_ss, comp)
}
comp_all_ss$model <- "single-species"

est_comp_ms <- comp_out(ms_estM1$model)
comp_all_ms <- data.frame()
for(i in 1:length(runs_ms)) {
  load(paste0("models/profile/ms/", runs_ms[i]))
  comp <- comp_out(run)
  comp$NLL <- comp$NLL - est_comp_ms$NLL
  comp$M1 <- round(run$quantities$M1[1, 1, 1], digits = 2)
  comp_all_ms <- rbind(comp_all_ms, comp)
}
comp_all_ms$model <- "cannibalism"

comp_all <- rbind(comp_all_ss, comp_all_ms)
comp_all$model <- factor(comp_all$model, levels = c("single-species", "cannibalism"))
comp_profile_plot <- ggplot(comp_all, aes(x = M1, y = NLL, color = component)) +
  geom_line() +
  geom_point(aes(shape = component)) +
  geom_vline(data = est_M1, mapping = aes(xintercept = value), linetype = "dashed") +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  ggsidekick::theme_sleek() +
  facet_wrap(~model)
comp_profile_plot

ggsave(filename="plots/CEATTLE/cannibalism/Testing/M1_comp_profile.png", 
       comp_profile_plot, 
       width=180, height=80, units="mm", dpi=300)

### Plot effect on spawning output --------------------------------------------
ssb_all_ss <- data.frame()
for(i in 1:length(runs_ss)) {
  load(paste0("models/profile/ss/", runs_ss[i]))
  ssb <- data.frame(t(run$quantities$biomassSSB))
  ssb$error <- (run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")] * 2)
  ssb$year <-rownames(ssb)
  rownames(ssb) <- NULL
  ssb$M1 <- round(run$quantities$M1[1, 1, 1], digits = 2)
  ssb_all_ss <- rbind(ssb_all_ss, ssb)
}
colnames(ssb_all_ss)[1] <- "SSB"
ssb_all_ss$model <- "single-species"

ssb_all_ms <- data.frame()
for(i in 1:length(runs_ms)) {
  load(paste0("models/profile/ms/", runs_ms[i]))
  ssb <- data.frame(t(run$quantities$biomassSSB))
  ssb$error <- (run$sdrep$sd[which(names(run$sdrep$value) == "biomassSSB")] * 2)
  ssb$year <-rownames(ssb)
  rownames(ssb) <- NULL
  ssb$M1 <- round(run$quantities$M1[1, 1, 1], digits = 2)
  ssb_all_ms <- rbind(ssb_all_ms, ssb)
}
colnames(ssb_all_ms)[1] <- "SSB"
ssb_all_ms$model <- "cannibalism"

ssb_all <- rbind(ssb_all_ss, ssb_all_ms)
ssb_all$SSB <- ssb_all$SSB / 1000000  # to Mt
ssb_all$error <- ssb_all$error / 1000000  # to Mt
ssb_all$year <- as.numeric(ssb_all$year)
ssb_all$model <- factor(ssb_all$model, levels = c("single-species", "cannibalism"))
ssb_all$M1 <- factor(as.character(ssb_all$M1))

ssb_profile_plot <- ggplot(ssb_all, aes(x = year, y = SSB, color = M1, fill = M1)) +
  geom_line() +
  geom_ribbon(aes(ymin=(SSB - error), ymax=(SSB + error)), alpha = 0.2, color = NA) + 
  ylab("Spawning Stock Biomass (Mt)") + labs(color = "M1") +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
  ggsidekick::theme_sleek() +
  facet_wrap(~model)
ssb_profile_plot

ggsave(filename="plots/CEATTLE/cannibalism/Testing/M1_profile_SSB.png", 
       ssb_profile_plot, 
       width=200, height=80, units="mm", dpi=300)
