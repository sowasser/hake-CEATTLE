# Profile over M1 

library(Rceattle)
library(ggplot2)

data <- read_data(file = "data/hake_intrasp_230427.xlsx")  # Read in data

### Run profile over M1 in single-species model -------------------------------
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
  save(run, file = paste0("models/profile/run", as.character(M1_change), ".Rdata"))
}

# Run the function for all values in M_vec
for(i in 1:length(M_vec)) {
  get_profile(M_vec[i], msm = 0)
}


# Update M_vec as needed if function/R crashes, or for more values
M_vec_new <- c(0.33)

# Run the function for all values in M_vec
for(i in 1:length(M_vec_new)) {
  get_profile(M_vec_new[i])
}

# Compare likelihoods & plot --------------------------------------------------
runs <- list.files(path = "models/profile")  # List of all model runs

# Get sum of joint negative log likelihood for each run
jnll_sums <- c() 
for(i in 1:length(runs)) {
  load(paste0("models/profile/", runs[i]))
  jnll_summary <- as.data.frame(run$quantities$jnll_comp)
  jnll_summary$sum <- rowSums(run$quantities$jnll_comp)
  jnll_sums[i] <- sum(rowSums(run$quantities$jnll_comp))
}

# Combine with M1 input for each run
profile <- cbind.data.frame(M1 = M_vec,
                            JNLL = jnll_sums)

# Convert to range from 0-1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

ggplot(profile, aes(x = M1, y = range01(JNLL))) +
  geom_line() +
  ylab("JNLL") +
  ggsidekick::theme_sleek()
