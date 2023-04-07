fit <- function(run) {
  objective <- run$opt$objective
  jnll <- run$quantities$jnll
  K <- run$opt$number_of_coefficients[1]
  AIC <- run$opt$AIC
  gradient <- run$opt$max_gradient
  
  fit <- cbind(objective, jnll, K, AIC, gradient)
  return(fit)
}


remove.packages("Rceattle")
remove.packages("00LOCK-Rceattle")
devtools::install_github("grantdadams/Rceattle")

HCR = Rceattle::build_hcr(HCR = 6, # Cat 1 HCR
                          FsprLimit = 0.4, # F40%
                          Ptarget = 0.4, # Target is 40% B0
                          Plimit = 0.1, # No fishing when SB<SB10
                          Pstar = 0.45,
                          Sigma = 0.5)

main_data <- Rceattle::read_data(file = "data/hake_intrasp_230324.xlsx")
main_data$est_M1 <- 0
main_run <- Rceattle::fit_mod(data_list = main_data,
                              msmMode = 0, # Single-species mode - no predation mortality
                              proj_mean_rec = 0,  # Project the model using: 0 = mean recruitment (average R of hindcast) or 1 = exp(ln_R0 + rec_devs)
                              estimateMode = 0,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
                              HCR = HCR,
                              phase = "default")
fit(main_run)
save(main_run, file = "main_run.Rdata")


remove.packages("Rceattle")
remove.packages("00LOCK-Rceattle")
devtools::install_github("grantdadams/Rceattle", ref = "dev")

dev_data <- Rceattle::read_data(file = "data/hake_intrasp_230324.xlsx")
dev_run <- Rceattle::fit_mod(data_list = dev_data,
                             msmMode = 0, # Single-species mode - no predation mortality
                             M1Fun = Rceattle::build_M1(M1_model = 0, 
                                                        updateM1 = TRUE),
                             estimateMode = 0,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
                             HCR = HCR,
                             phase = "default",
                             initMode = 1)

fit(dev_run)

save(dev_run, file = "dev_run.Rdata")

load("main_run.Rdata")

inits <- main_run$estimated_params
inits$rec_pars <- matrix(9, nrow = main_run$data_list$nspp, ncol = 2)
inits$rec_pars[,1] <- inits$ln_mean_rec

inits$ln_mean_rec <- NULL

##########################################
ss_run <- Rceattle::fit_mod(data_list = dev_data,
                            file = NULL,
                            inits = inits, # Initial parameters = 0
                            estimateMode = 3, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            verbose = 1,
                            phase = "default",
                            M1Fun = Rceattle::build_M1(M1_model = 1, 
                                                       updateM1 = FALSE),
                            initMode = 1)

ss_run$quantities$jnll_comp
