# Run CEATTLE for hake data with no diet data 
devtools::install_github("grantdadams/Rceattle@dev")
library(Rceattle)

library(r4ss)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(viridis)

# Run CEATTLE -----------------------------------------------------------------
# mydata <- Rceattle::read_data( file = "data/hake_from_ss.xlsx")
hake_nodiet <- Rceattle::read_data(file = "data/hake_singlesp_220401.xlsx")

nodiet_run <- Rceattle::fit_mod(data_list = hake_nodiet,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                            # debug = 1, # 1 = estimate, 0 = don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              phase = "default")

# Check what all comes out of CEATTLE
ceattle_stuff <- nodiet_run$quantities

years <- 1980:2022

# Pull out SSB & total biomass from CEATTLE & combine -------------------------
nodiet_ssb <- (c(nodiet_run$quantities$biomassSSB) * 2)
nodiet_biomass <- c(nodiet_run$quantities$biomass)

nodiet_biom_wide <- as.data.frame(cbind(years, nodiet_ssb, nodiet_biomass))
colnames(nodiet_biom_wide) <- c("year", "SSB", "Total Biomass")
nodiet_biom <- melt(nodiet_biom_wide, id.vars = "year")
colnames(nodiet_biom)[2:3] <- c("type", "CEATTLE no diet")

write.csv(nodiet_biom, "data/ceattle_nodiet_biom.csv", row.names = FALSE)


# Recruitment -----------------------------------------------------------------
nodiet_R <- c(nodiet_run$quantities$R)
write.csv(nodiet_R, "data/ceattle_nodiet_R.csv", row.names = FALSE)


# Numbers at age --------------------------------------------------------------
nbyage_raw <- as.data.frame(as.table(nodiet_run$quantities$NByage))

# Remove rows with 0s, artifact of a single-sex model
nbyage <- nbyage_raw[-seq(0, nrow(nbyage_raw), 2), -c(1:2)]

levels(nbyage$Var3) <- c(1:20)
levels(nbyage$Var4) <- c(1980:2022)
colnames(nbyage) <- c("age", "year", "numbers")

write.csv(nbyage, "data/ceattle_nodiet_nbyage.csv", row.names = FALSE)

# Read in and transform data from stock synthesis - with & without age0 and for
# beginning and middle of the year to find best match for CEATTLE
nbyage_ss_all <- read.csv("data/assessment/nbyage.csv")

b_wide <- nbyage_ss_all[1:43, -c(2:3)]
m_wide <- nbyage_ss_all[44:86, -c(2:3)]

new_nbyage <- function(df, name) {
  colnames(df) <- c("year", as.character(1:20))
  df2 <- melt(df, id.vars = "year")
  df3 <- cbind(df2, rep(name, nrow(df2)))
  colnames(df3)[2:4] <- c("age", "numbers", "source")
  
  return(df3)
}

b <- new_nbyage(b_wide, "SS - beginning")
m <- new_nbyage(m_wide, "SS - middle")

nbyage_ss <- rbind(b, m)

# Combine with CEATTLE & plot both per-year and mean
nbyage_ss <- nbyage_ss[, c(2, 1, 3, 4)]
nbyage_ceattle <- cbind(nbyage, rep("CEATTLE", nrow(nbyage)))
colnames(nbyage_ceattle)[4] <- "source"

nbyage_final <- rbind(nbyage_ceattle, nbyage_ss)

nbyage_plot_year <- ggplot(nbyage_final, aes(x=age, y=numbers, fill=source)) +
  geom_bar(stat = "identity", position = "dodge", ) +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("age") + ylab("numbers") +
  facet_wrap(~year)
nbyage_plot_year

ggsave(filename = "plots/CEATTLE/nbyage_year.png", nbyage_plot_year, 
       width=500, height=400, units="mm", dpi=300)

nbyage_mean <- nbyage_final %>% group_by(age, source) %>%
  summarise(mean_number = mean(numbers))

nbyage_plot_mean <- ggplot(nbyage_mean, aes(x=age, y=mean_number, fill=source)) +
  geom_bar(stat = "identity", position = "dodge", ) +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  xlab("age") + ylab("numbers") 
nbyage_plot_mean

ggsave(filename = "plots/CEATTLE/nbyage_mean.png", nbyage_plot_mean, 
       width=200, height=120, units="mm", dpi=300)


# # Run r4ss diagnostics on stock synthesis model -------------------------------
# # Diagnostics from r4ss
# mydir <- file.path(file.path("hake_assessment/2020_Hake_Assessment"))
# 
# # read the model output and print diagnostic messages 
# replist <- SS_output(dir = mydir, 
#                      verbose = TRUE,
#                      printstats = TRUE)
# 
# # plots the results
# SS_plots(replist)
