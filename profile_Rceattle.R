# Profile over M1 

library(Rceattle)
library(ggplot2)
library(dplyr)

data <- read_data(file = "data/hake_intrasp_230427.xlsx")  # Read in data

M_vec <- seq(0.15, 0.35, 0.01)  # Initial vector of M values to test

# Function updating M1 for each run of the model
get_profile <- function(M1_change, msm) {
  data$M1_base[, 3:17] <- M1_change
  run <- fit_mod(
    data_list = data,
    inits = NULL,
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

### Run profile over M1 in single-species model -------------------------------
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

# Get JNLL for each run 
runs <- list.files(path = "models/profile/ss")  # List of all model runs

# Get sum of joint negative log likelihood for each run
jnll_sums_ss <- c() 
for(i in 1:length(runs)) {
  load(paste0("models/profile/ss/", runs[i]))
  jnll_summary <- as.data.frame(run$quantities$jnll_comp)
  jnll_summary$sum <- rowSums(run$quantities$jnll_comp)
  jnll_sums_ss[i] <- sum(rowSums(run$quantities$jnll_comp))
}

# Combine with M1 input for each run
profile_ss <- cbind.data.frame(M1 = M_vec,
                               JNLL = jnll_sums_ss,
                               model = "single-species")


### Run profile over M1 in cannibalism model ----------------------------------
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

# Get JNLL for each run & plot with M1 value
runs <- list.files(path = "models/profile/ms")  # List of all model runs

# Get sum of joint negative log likelihood for each run
jnll_sums_ms <- c() 
for(i in 1:length(runs)) {
  load(paste0("models/profile/ms/", runs[i]))
  jnll_summary <- as.data.frame(run$quantities$jnll_comp)
  jnll_summary$sum <- rowSums(run$quantities$jnll_comp)
  jnll_sums_ms[i] <- sum(rowSums(run$quantities$jnll_comp))
}

# Combine with M1 input for each run
profile_ms <- cbind.data.frame(M1 = M_vec,
                               JNLL = jnll_sums_ms,
                               model = "cannibalism")


### Plot profile --------------------------------------------------------------
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

ggsave(filename="plots/CEATTLE/cannibalism/M1_profile.png", 
       profile_plot, 
       width=180, height=80, units="mm", dpi=300)
