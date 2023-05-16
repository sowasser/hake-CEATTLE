# Profile over M1 

library(Rceattle)
library(ggplot2)
library(viridis)
library(dplyr)

data <- read_data(file = "data/hake_intrasp_230427.xlsx")  # Read in data
load("models/ss_estM1.Rdata")
load("models/ms_estM1.Rdata")

startM_ss <- round(ss_estM1$model$quantities$M1[1, 1, 1], digits = 2)
startM_ms <- round(ms_estM1$model$quantities$M1[1, 1, 1], digits = 2)
nll_ss <- ss_estM1$model$quantities$jnll
nll_ms <- ms_estM1$model$quantities$jnll

delta <- seq(0.01,0.10, 0.01)

M_ss <- rep(round(ss_estM1$model$quantities$M1[1, 1, 1], digits = 2), 10)

M_vec_ss <- c(c(M_ss - delta), 
              round(ss_estM1$model$quantities$M1[1, 1, 1], digits = 2),
              c(M_ss + delta))



# Function updating M1 for each run of the model
get_profile <- function(M1_change, init, msm) {
  data$M1_base[, 3:17] <- M1_change
  run <- fit_mod(
    data_list = data,
    inits = init,
    file = NULL, # Don't save
    msmMode = msm, 
    M1Fun = Rceattle::build_M1(M1_model = 0,
                               updateM1 = FALSE,
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
}

### Run profile over M1 -------------------------------------------------------
# SINGLE SPECIES 
# Run the function for all values in M_vec
for(i in 1:length(M_vec)) {
  get_profile(M_vec[i], msm = 0)
}

# Update M_vec as needed if function/R crashes, or for more values
M_vec_new <- c(0.33)

# Run the function for all values in M_vec
for(i in 1:length(M_vec_new)) {
  get_profile(M_vec_new[i], msm = 0)
}

# CANNIBALISM
# Run the function for all values in M_vec
for(i in 1:length(M_vec)) {
  get_profile(M_vec[i], msm = 1)
}

# Update M_vec as needed if function/R crashes, or for more values
M_vec_new <- c(0.24)

# Run the function for all values in M_vec
for(i in 1:length(M_vec_new)) {
  get_profile(M_vec_new[i], msm = 1)
}


### Get JNLL for each run and plot --------------------------------------------
M_vec <- M_vec <- seq(0.15, 0.35, 0.01) 

# Single species
runs_ss <- list.files(path = "models/profile/ss")  # List of all model runs

# Get sum of joint negative log likelihood for each run
jnll_all_ss <- c() 
for(i in 1:length(runs_ss)) {
  load(paste0("models/profile/ss/", runs_ss[i]))
  jnll_all_ss[i] <- run$quantities$jnll
}

# Combine with M1 input for each run
profile_ss <- cbind.data.frame(M1 = M_vec,
                               JNLL = jnll_all_ss - nll_ss,
                               model = "single-species")

# Cannibalism
runs_ms <- list.files(path = "models/profile/ms")  # List of all model runs

# Get sum of joint negative log likelihood for each run
jnll_all_ms <- c() 
for(i in 1:length(runs_ms)) {
  load(paste0("models/profile/ms/", runs_ms[i]))
  jnll_all_ms[i] <- run$quantities$jnll
}

# Combine with M1 input for each run
profile_ms <- cbind.data.frame(M1 = M_vec,
                               JNLL = jnll_all_ms - nll_ms,
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
  ggsidekick::theme_sleek() +
  facet_wrap(~model)
profile_plot

ggsave(filename="plots/CEATTLE/cannibalism/Testing/M1_profile.png", 
       profile_plot, 
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
