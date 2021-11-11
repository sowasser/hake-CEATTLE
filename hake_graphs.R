# Exploring the hake dynamics from the stock assessment data:
# https://github.com/pacific-hake/hake-assessment/tree/master/data

library(ggplot2)
library(ggsidekick)
library(reshape2)
library(viridis)
library(dplyr)

path <- "~/Desktop/Local/hake-assessment-master/data/"

# Weight & length at age ------------------------------------------------------
maturity_table <- read.csv(paste0(path, "maturity-table.csv"))
hake_maturity_data <- read.csv(paste0(path, "hake-maturity-data.csv"))

# Weight at age
weight_age <- ggplot(maturity_table, aes(x=age, y=avg.wt)) +
  geom_point() +
  theme_sleek()
# weight_age

ggsave(filename="plots/weight_age.pdf", weight_age,
       width=150, height=100, units="mm", dpi=300)

# Yearly weight at age
weight_age_yearly <- ggplot(hake_maturity_data, aes(x=Age, y=Weight_kg)) +
  geom_point() +
  theme_sleek() +
  facet_wrap(~Year, ncol = 3) 
# weight_age_yearly

ggsave(filename="plots/weight_age_yearly.pdf", weight_age_yearly,
       width=300, height=200, units="mm", dpi=300)

# Yearly length at age
length_age_yearly <- ggplot(hake_maturity_data, aes(x=Age, y=Length_cm)) +
  geom_point() +
  theme_sleek() +
  facet_wrap(~Year, ncol = 3) 
# length_age_yearly

ggsave(filename="plots/length_age_yearly.pdf", length_age_yearly,
       width=300, height=200, units="mm", dpi=300)

# Weight/length relationship (sanity check)
weight_length <- ggplot(hake_maturity_data, aes(x=Length_cm, y=Weight_kg)) +
  geom_point() +
  theme_sleek()
# weight_length

# ggsave(filename="plots/weight_length.pdf", weight_length,
#        width=150, height=100, units="mm", dpi=300)


# Maturity --------------------------------------------------------------------
# Age and maturity
age_maturity <- ggplot(maturity_table, aes(x=age, y=maturity)) +
  geom_point() +
  theme_sleek() 
# age_maturity

ggsave(filename="plots/age_maturity.pdf", age_maturity,
       width=150, height=100, units="mm", dpi=300)

# Maturity ogives
maturity_ogives_wide <- read.csv(paste0(path, "maturity-ogives.csv"))
maturity_ogives <- melt(maturity_ogives_wide, id.vars = "length.cm",
                        variable.name = "source", value.name = "ogive")

ogives <- ggplot(maturity_ogives, aes(x=length.cm, y=ogive, color=source)) +
  geom_point() +
  theme_sleek() 
# ogives


# Age composition -------------------------------------------------------------
us_cp_age_data <- read.csv(paste0(path, "us-cp-age-data.csv"))
us_ms_age_data <- read.csv(paste0(path, "us-ms-age-data.csv"))
us_shore_age_data <- read.csv(paste0(path, "us-shore-age-data.csv"))

us_cp_age_data <- cbind(us_cp_age_data[, -3], source=rep("cp", length(us_cp_age_data[, 1])))
us_ms_age_data <- cbind(us_ms_age_data[, -3], source=rep("ms", length(us_ms_age_data[, 1])))
us_shore_age_data <- cbind(us_shore_age_data[, -3], source=rep("shore", length(us_shore_age_data[, 1])))
age_data_wide <- rbind(us_cp_age_data, us_ms_age_data, us_shore_age_data)
age_data <- melt(age_data_wide, id.vars = c("year", "n.fish", "source"),
                 variable.name = "age", value.name = "proportion")
age_data <- cbind(age_data, catch=(age_data$n.fish * age_data$proportion))

catch_age_comp <- ggplot(age_data, aes(x=year, y=catch, fill=age)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~source)
# catch_age_comp

ggsave(filename="plots/catch_age_comp.pdf", catch_age_comp,
       width=300, height=100, units="mm", dpi=300)


# Catch rates -----------------------------------------------------------------
# Import data, rename columns, add column with country & source names, reshape
catch_data <- function(file, country, type) {
  df <- read.csv(paste0(path, file))
  
  if (country == "Canada") {
    colnames(df) <- c("year", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
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
us_unid <- as.data.frame(cbind(year = us_unid[, 3], nation = rep("US", length(us_unid[, 3])), 
                               source = rep("unidentified", length(us_unid[, 3])), month = us_unid[, 2], 
                               catch = us_unid[, 4]))

monthly_catch <- rbind(catch_data("can-ft-catch-by-month.csv", "Canada", "freezer-trawler"),
                       catch_data("can-jv-catch-by-month.csv", "Canada", "JV"),
                       catch_data("can-ss-catch-by-month.csv", "Canada", "shoreside"),
                       catch_data("us-cp-catch-by-month.csv", "US", "catcher-processor"),
                       catch_data("us-ms-catch-by-month.csv", "US", "mothership"),
                       catch_data("us-research-catch-by-month.csv", "US", "research"),
                       catch_data("us-shore-catch-by-month.csv", "US", "shore-based"),
                       us_unid)
monthly_catch$catch <- as.numeric(monthly_catch$catch)
all_catch <- monthly_catch %>% group_by(year, nation, source) %>%
  summarize(total_catch = sum(catch))

catch <- ggplot(all_catch, aes(y=total_catch, x=year, fill=source)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~nation, ncol=2)

ggsave(filename="plots/catch.pdf", catch,
       width=400, height=200, units="mm", dpi=300)


# Survey ----------------------------------------------------------------------
survey_history <- read.csv(paste0(path, "survey-history.csv"))
survey_history$year <- as.factor(survey_history$year)

survey_biomass <- ggplot(survey_history, aes(y=biomass, x=year)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=biomass-cv, ymax=biomass+cv), width=.2) +
  theme_sleek()
# survey_biomass

ggsave(filename="plots/survey_biomass.pdf", survey_biomass,
       width=150, height=100, units="mm", dpi=300)
