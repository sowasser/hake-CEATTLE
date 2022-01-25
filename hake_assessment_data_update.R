# Script for data manipulation for adapting the hake stock assessment for 
# input into CEATTLE

library(reshape2)
library(dplyr)

path <- "~/Desktop/Local/hake-assessment-master/data/"

# Catch data ------------------------------------------------------------------
# Comes from multiple sources in the US and Canada

# Import data, rename columns, add column with country & source names, reshape
catch_data <- function(file, country, type) {
  df <- read.csv(paste0(path, file))
  
  if (country == "Canada") {
    colnames(df) <- c("year", 1:12)
    df2 <- cbind(df, nation = rep("Canada", length(df[, 1])), 
                 source = rep(type, length(df[, 1])))
    df3 <- melt(df2, id.vars = c("year", "nation", "source"), 
                variable.name = "month", value.name = "catch")
    return(df3)
  }
  
  if (country == "US") {
    df2 <- cbind(year = df[, 2], nation = rep("US", length(df[, 2])), 
                 source = rep(type, length(df[, 2])), month = df[, 1], 
                 catch = df[, 3])
    return(as.data.frame(df2))
  }
  
}

us_unid <- read.csv(paste0(path, "us-catch-by-month.csv"))
us_unid <- as.data.frame(cbind(year = us_unid[, 3], 
                               nation = rep("US", length(us_unid[, 3])), 
                               source = rep("unidentified", length(us_unid[, 3])),
                               month = us_unid[, 2], 
                               catch = us_unid[, 4]))

catch_all <- rbind(catch_data("can-ft-catch-by-month.csv", "Canada", "freezer-trawler"),
                   catch_data("can-jv-catch-by-month.csv", "Canada", "JV"),
                   catch_data("can-ss-catch-by-month.csv", "Canada", "shoreside"),
                   catch_data("us-cp-catch-by-month.csv", "US", "catcher-processor"),
                   catch_data("us-ms-catch-by-month.csv", "US", "mothership"),
                   catch_data("us-research-catch-by-month.csv", "US", "research"),
                   catch_data("us-shore-catch-by-month.csv", "US", "shoreside"))
                 # us_unid)

catch_all$catch <- as.numeric(catch_all$catch)

monthly_catch <- catch_all %>% group_by(nation, year, month) %>%
  summarize(total_catch = sum(catch))

write.csv(monthly_catch, "data/assessment/monthly_catch.csv")

catch_source <- catch_all %>% group_by(year, nation, source) %>%
  summarize(total_catch = sum(catch))

write.csv(catch_all, "data/assessment/catch_source.csv")


# Age composition -------------------------------------------------------------
# US data
us_cp_age <- read.csv(paste0(path, "us-cp-age-data.csv"))
us_ms_age <- read.csv(paste0(path, "us-ms-age-data.csv"))
us_shore_age <- read.csv(paste0(path, "us-shore-age-data.csv"))

us_cp_age <- cbind(us_cp_age, source=rep("cp", length(us_cp_age[, 1])))
us_ms_age <- cbind(us_ms_age, source=rep("ms", length(us_ms_age[, 1])))
us_shore_age <- cbind(us_shore_age, source=rep("shore", length(us_shore_age[, 1])))
colnames(us_shore_age)[3] <- "n.hauls"

us_age_wide <- rbind(us_cp_age, us_ms_age, us_shore_age)
us_age <- melt(us_age_wide, id.vars = c("year", "n.fish", "n.hauls", "source"),
               variable.name = "age", value.name = "proportion")

# # Canadian data
# can_age_all <- read.csv(paste0(path, "can-age-data.csv"), skip=1)
# can_shore_age <- cbind(can_age_all[1:27,], can_age_all[65:91,2], rep("shore", length(27)))
# colnames(can_shore_age) <- c("year", 1:15, "n.hauls", "source")
# can_ft_age <- cbind(can_age_all[29:43,], can_age_all[93:107, 2], rep("ft", length(15)))
# colnames(can_ft_age) <- c("year", 1:15, "n.hauls", "source")
# can_jv_age <- cbind(can_age_all[45:63,], can_age_all[109:127, 2], rep("jv", length(19)))
# colnames(can_jv_age) <- c("year", 1:15, "n.hauls", "source")
# 
# can_age_wide <- rbind(can_shore_age, can_ft_age, can_jv_age)
# can_age <- melt(can_age_wide, id.vars = c("year", "n.hauls", "source"),
#                 variable.name = "age", value.name = "proportion")

us_age <- cbind(us_age, catch=(us_age$n.fish * us_age$proportion))

write.csv(us_age, "data/assessment/age_comp.csv")


# Get wide format age composition for NByageFixed - filter age comp by year, 
# find mean proportion, multiply by number of fish
n_year <- list()
for(i in 2008:2020) {
  df <- us_age_wide %>% filter(year == i)
  avg <- sapply(df[4:18], mean)
  number <- avg * df$n.fish
  n_year[[i]] <- number
}
n_final <- do.call(rbind, n_year)

n_final <- as.data.frame(cbind(2008:2020, n_final))
colnames(n_final)[1] <- "year"

write.csv(n_final, "data/assessment/age_comp_yearly.csv")


# Weight/length at age & maturity ---------------------------------------------
# Attempt to get wide format for weight at age - will only be empirical
hake_maturity_data <- read.csv(paste0(path, "hake-maturity-data.csv"))
weight_age_long <- na.omit(hake_maturity_data[, c(2, 10, 12)])

weight_age_wide <- dcast(weight_age_long, Year ~ Age, value.var = "Weight_kg", 
                         fun.aggregate = mean)

write.csv(weight_age_wide, "data/assessment/empirical_weight_age.csv")

# Length at age
length_age_all <- na.omit(hake_maturity_data[, 11:12])
lengths <- unique(length_age_all$Length_cm)
mean_lengths <- length_age_all %>% group_by(Age) %>% 
  summarise(Length = mean(Length_cm))

write.csv(mean_lengths, "data/assessment/length_at_age.csv")
