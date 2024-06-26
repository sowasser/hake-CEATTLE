# Profile over M1 

# Re-install Rceattle if needed
# remove.packages("Rceattle")
# remove.packages("00LOCK-Rceattle")
# devtools::install_github("grantdadams/Rceattle", ref = "dev")

library(Rceattle)
library(ggplot2)
library(viridis)
library(dplyr)
# Set ggplot theme
theme_set(ggsidekick::theme_sleek())

data <- read_data(file = "data/hake_intrasp_240102.xlsx")  # Read in data

# Function updating M1 for each run of the model
get_profile <- function(new_M1, model, msm) {
  # M1_change <- -0.01
  # model <- ss_estM1$model$estimated_params
  # msm <- 0
  M1_base <- data$M1_base[, 3:22]
  data$M1_base[, 3:22] <- new_M1
  run <- fit_mod(
    data_list = data,
    inits = model,
    msmMode = msm,
    M1Fun = Rceattle::build_M1(M1_model = 0,
                               updateM1 = TRUE,
                               M1_use_prior = FALSE),
    estimateMode = 1,  # 0 = Fit the hindcast model and projection with HCR specified via HCR
    phase = "default",
    initMode = 1,
    loopnum = 7
  )

  # Save the resulting run to losing progress to R bombs!
  if (msm == 0) {save(run, file = paste0("models/profile/ss/run",
                                         as.character(new_M1),
                                         ".Rdata"))}
  if (msm == 1) {save(run, file = paste0("models/profile/ms/run",
                                         as.character(new_M1),
                                         ".Rdata"))}
  
  message(run$opt$obj)  # print objective function to console
  message(run$quantities$jnll)  # print JNLL to console
  return(run)
}

### Run profile over M1 -------------------------------------------------------
# Load model with estimated M1 & check starting value
load("models/ss_estM1.Rdata")
startM_ss <- round(exp(ss_estM1$model$initial_params$ln_M1)[1, 1, 1], digits = 2)

run0 <- get_profile(startM_ss, ss_estM1$model$estimated_params, 0)  # 0.24
# SS down
run1 <- get_profile(0.23, run0$estimated_params, 0)
run2 <- get_profile(0.22, run1$estimated_params, 0)
run3 <- get_profile(0.21, run2$estimated_params, 0)
run4 <- get_profile(0.20, run3$estimated_params, 0)
run5 <- get_profile(0.19, run4$estimated_params, 0)
run6 <- get_profile(0.18, run5$estimated_params, 0)
run7 <- get_profile(0.17, run6$estimated_params, 0)
run8 <- get_profile(0.16, run7$estimated_params, 0)
run9 <- get_profile(0.15, run8$estimated_params, 0)

# SS up
run10 <- get_profile(0.25, run0$estimated_params, 0)
run11 <- get_profile(0.26, run10$estimated_params, 0)
run12 <- get_profile(0.27, run11$estimated_params, 0)
run13 <- get_profile(0.28, run12$estimated_params, 0)
run14 <- get_profile(0.29, run13$estimated_params, 0)
run15 <- get_profile(0.30, run14$estimated_params, 0)

rm(list = ls())  # clear environment to re-set runs
# Load model with estimated M1 & check starting value
load("models/ms_estM1.Rdata")
startM_ms <- round(exp(ms_estM1$model$initial_params$ln_M1)[1, 1, 1], digits = 2)

run0 <- get_profile(startM_ms, ms_estM1$model$estimated_params, 1)  # 0.27

# ms down
run1 <- get_profile(0.26, run0$estimated_params, 1)
run2 <- get_profile(0.25, run1$estimated_params, 1)
run3 <- get_profile(0.24, run2$estimated_params, 1)
run4 <- get_profile(0.23, run3$estimated_params, 1)
run5 <- get_profile(0.22, run4$estimated_params, 1)
run6 <- get_profile(0.21, run5$estimated_params, 1)
run7 <- get_profile(0.20, run6$estimated_params, 1)
run8 <- get_profile(0.19, run7$estimated_params, 1)
run9 <- get_profile(0.18, run8$estimated_params, 1)
run10 <- get_profile(0.17, run9$estimated_params, 1)
run11 <- get_profile(0.16, run10$estimated_params, 1)
run12 <- get_profile(0.15, run11$estimated_params, 1)

# ms up
run13 <- get_profile(0.28, run0$estimated_params, 1)
run14 <- get_profile(0.29, run13$estimated_params, 1)
run15 <- get_profile(0.30, run14$estimated_params, 1)


### Get JNLL for each run and plot --------------------------------------------
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
                               JNLL = jnll_all_ss,
                               relative_NLL = jnll_all_ss - min(jnll_all_ss),
                               model = "single-species") %>%
  arrange(M1)

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
profile_ms <- cbind.data.frame(M1 = m1_all_ss,
                               JNLL = jnll_all_ms,
                               relative_NLL = jnll_all_ms - min(jnll_all_ms),
                               model = "cannibalism")


# add points at estimated value
load("models/ss_estM1.Rdata")
load("models/ms_estM1.Rdata")
est_points <- cbind.data.frame(model = factor(c("single-species", "cannibalism")),
                               M1 = c(round(ss_estM1$model$quantities$M1[1, 1, 1], digits = 2),
                                      round(ms_estM1$model$quantities$M1[1, 1, 1], digits = 2)),
                               JNLL = c(ss_estM1$model$quantities$jnll - min(jnll_all_ss), 
                                        ms_estM1$model$quantities$jnll - min(jnll_all_ms)))

# add points at prior value
load("models/ss_priorM1.Rdata")
load("models/ms_priorM1.Rdata")
prior_points <- cbind.data.frame(model = factor(c("single-species", "cannibalism")),
                                 M1 = c(round(ss_priorM1$model$quantities$M1[1, 1, 1], digits = 2),
                                        round(ms_priorM1$model$quantities$M1[1, 1, 1], digits = 2)),
                                 JNLL = c(ss_priorM1$model$quantities$jnll - min(jnll_all_ss), 
                                          ms_priorM1$model$quantities$jnll - min(jnll_all_ms)))

all_profile <- rbind(profile_ss, profile_ms) %>%
  mutate(model = factor(model, levels = c("single-species", "cannibalism")))

# profile_plot <- ggplot() +
#   geom_line(data = all_profile, aes(x = M1, y = relative_NLL), linewidth = 1) +
#   geom_point(data = est_points, aes(x = M1, y = JNLL, color = model), size = 5, shape = "circle") +
#   geom_point(data = prior_points, aes(x = M1, y = JNLL, color = model), size = 5, shape = "triangle") +
#   geom_hline(yintercept = 2, color = "lightgray") +
#   # geom_vline(data = est_M1, mapping = aes(xintercept = value, color = model),
#   #            linetype = "dashed", linewidth = 1) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.5) +
#   ylab("JNLL") + labs(color = "M1 estimate") +
#   # scale_x_continuous(breaks = seq(0.15, 0.35, 0.01)) +  # all scale markers for investigating
#   ggsidekick::theme_sleek() +
#   facet_wrap(~model)
# profile_plot

# ggsave(filename="plots/CEATTLE/cannibalism/Testing/M1/M1_profile.png",
#        profile_plot,
#        width=180, height=80, units="mm", dpi=300)

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
  # # Select components for easier plotting.
  comp <- comp %>% filter(!component %in% c("Fishery age composition", 
                                            "Survey age composition", 
                                            "Total catch"))
  return(comp)
}

load("models/profile/ss/run0.22.Rdata"); ss_low <- run
comp_all_ss <- data.frame()
for(i in 1:length(runs_ss)) {
  load(paste0("models/profile/ss/", runs_ss[i]))
  comp <- comp_out(run)
  comp$NLL <- comp$NLL - comp_out(ss_low)$NLL
  comp$M1 <- round(run$quantities$M1[1, 1, 1], digits = 2)
  comp_all_ss <- rbind(comp_all_ss, comp)
}
comp_all_ss$model <- "Single-Species"

load("models/profile/ms/run0.26.Rdata"); ms_low <- run
comp_all_ms <- data.frame()
comp_all_ms <- data.frame()
for(i in 1:length(runs_ms)) {
  load(paste0("models/profile/ms/", runs_ms[i]))
  comp <- comp_out(run)
  comp$NLL <- comp$NLL - comp_out(ms_low)$NLL
  comp$M1 <- round(run$quantities$M1[1, 1, 1], digits = 2)
  comp_all_ms <- rbind(comp_all_ms, comp)
}
comp_all_ms$model <- "Cannibalism"

comp_all <- rbind(comp_all_ss, comp_all_ms)
comp_all$model <- factor(comp_all$model, levels = c("Single-Species", "Cannibalism"))

comp_all$component[comp_all$component == "Survey biomass"] <- "Survey"
comp_all$component[comp_all$component == "Age/length composition data"] <- "Age"
comp_all$component[comp_all$component == "Selectivity deviates"] <- "Sel"
comp_all$component[comp_all$component == "Recruitment deviates"] <- "Rec"
comp_all$component[comp_all$component == "Initial abundance deviates"] <- "Init"
comp_all$component[comp_all$component == "Fishing mortality deviates"] <- "F"

comp_profile_plot <- ggplot() +
  geom_line(data = (comp_all %>% filter(component == "Total NLL")), 
            aes(x = M1, y = NLL, color = component), linewidth = 1) +
  geom_line(data = (comp_all %>% filter(component != "Total NLL")), 
            aes(x = M1, y = NLL, color = component)) +
  geom_point(data = (comp_all), 
             aes(x = M1, y = NLL, shape = component, color = component), size = 2) +
  geom_hline(yintercept = 2, color = "lightgray") +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  ylab("Change in Negative Log Likelihood") +
  labs(shape = "Component", color = "Component") +
  ylim(-13.7, 13.7) +
  facet_wrap(~model, ncol = 1)
comp_profile_plot

ggsave(filename="~/Desktop/M1_comp_profile.tiff",
       comp_profile_plot,
       width=85, height=100, units="mm", dpi=500)

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
ssb_all_ss$model <- "Single-Species"

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
ssb_all_ms$model <- "Cannibalism"

ssb_all <- rbind(ssb_all_ss, ssb_all_ms)
ssb_all$SSB <- ssb_all$SSB / 1000000  # to Mt
ssb_all$error <- ssb_all$error / 1000000  # to Mt
ssb_all$year <- as.numeric(ssb_all$year)
ssb_all$model <- factor(ssb_all$model, levels = c("Single-Species", "Cannibalism"))
ssb_all$M1 <- factor(as.character(ssb_all$M1))

ssb_all_est <- ssb_all %>% filter(model == "Single-Species" & M1 == 0.21 | 
                                    model == "Cannibalism" & M1 == 0.23)

ssb_profile_plot <- ggplot() +
  geom_line(data = ssb_all, aes(x = year, y = SSB, color = M1)) +
  geom_ribbon(data = ssb_all, 
              aes(x = year, y = SSB, ymin=(SSB - error), ymax=(SSB + error), fill = M1), 
              alpha = 0.2, color = NA) + 
  geom_line(data = ssb_all_est, aes(x = year, y = SSB), linewidth = 1) +
  geom_vline(xintercept = 2020, linetype = 2, colour = "gray") +  # Add line at end of hindcast
  xlim(1980, 2022) +
  ylab("Spawning Biomass (Mt)") + xlab("Year") + labs(color = "M1") +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
  ggsidekick::theme_sleek() +
  facet_wrap(~model)
ssb_profile_plot

ggsave(filename="plots/CEATTLE/cannibalism/Testing/M1/M1_profile_SSB.png", 
       ssb_profile_plot, 
       width=200, height=80, units="mm", dpi=300)

### Plot fit to survey index --------------------------------------------------
init_surv <- run$data_list$srv_biom %>% filter(Year > 1)  # input survey biomass

survey_biom <- function(run, name) {
  srv <- data.frame(year = 1995:2019,
                    biomass = run$quantities$srv_bio_hat,
                    log_sd = run$quantities$srv_log_sd_hat,
                    M1 = rep(name, length(1995:2019)))
  return(srv)
}

ss_srv <- data.frame()
for(i in 1:length(runs_ss)) {
  load(paste0("models/profile/ss/", runs_ss[i]))
  srv_out <- survey_biom(run = run, 
                         name = round(run$quantities$M1[1, 1, 1], digits = 2))
  ss_srv <- rbind(ss_srv, srv_out)
}
ss_srv$model <- "single-species"

ms_srv <- data.frame()
for(i in 1:length(runs_ms)) {
  load(paste0("models/profile/ms/", runs_ms[i]))
  srv_out <- survey_biom(run = run, 
                         name = round(run$quantities$M1[1, 1, 1], digits = 2))
  ms_srv <- rbind(ms_srv, srv_out)
}
ms_srv$model <- "cannibalism"

assess_srv <- read.csv(paste0("data/assessment/2020/survey_out.csv"))

srv_all <- rbind(ss_srv, ms_srv)
srv_all$model <- factor(srv_all$model, levels = c("single-species", "cannibalism"))
srv_all$M1 <- factor(as.character(srv_all$M1))

survey_profile_plot <- ggplot() +
  geom_pointrange(data = init_surv, 
                  aes(x = Year, y = Observation,
                      ymin = exp(log(Observation) - 1.96*Log_sd),
                      ymax = exp(log(Observation) + 1.96*Log_sd)),
                  fatten = 5) +
  geom_line(data = srv_all, aes(x = year, y = biomass, color = M1), alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +
  geom_line(data = assess_srv, aes(x = Yr, y = Exp), linetype = "dashed") +
  xlab("year") + ylab("Index of Abundance") +
  ylim(0, NA) +
  facet_wrap(~model, ncol = 2)
survey_profile_plot

ggsave(filename="plots/CEATTLE/cannibalism/Testing/M1/M1_profile_survey.png", 
       survey_profile_plot, 
       width=200, height=80, units="mm", dpi=300)
