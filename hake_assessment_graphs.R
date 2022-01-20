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

# Overall weight at age
weight_age <- ggplot(maturity_table, aes(x=age, y=avg.wt)) +
  geom_point() +
  theme_sleek()
# weight_age

ggsave(filename="plots/assessment/weight_age.pdf", weight_age,
       width=150, height=100, units="mm", dpi=300)

# Yearly weight at age
weight_age_yearly <- ggplot(hake_maturity_data, aes(x=Age, y=Weight_kg)) +
  geom_point() +
  theme_sleek() +
  facet_wrap(~Year, ncol = 3) 
# weight_age_yearly

ggsave(filename="plots/assessment/weight_age_yearly.pdf", weight_age_yearly,
       width=300, height=200, units="mm", dpi=300)

# Overall length at age
length_age <- ggplot(hake_maturity_data, aes(x=Length_cm, y=Age)) +
  geom_point() +
  theme_sleek() 
# length_age

ggsave(filename="plots/assessment/length_age.pdf", length_age,
       width=150, height=100, units="mm", dpi=300)

# Yearly length at age
length_age_yearly <- ggplot(hake_maturity_data, aes(x=Length_cm, y=Age)) +
  geom_point() +
  theme_sleek() +
  facet_wrap(~Year, ncol = 3) 
# length_age_yearly

ggsave(filename="plots/assessment/length_age_yearly.pdf", length_age_yearly,
       width=300, height=200, units="mm", dpi=300)

# Weight/length relationship (sanity check)
weight_length <- ggplot(hake_maturity_data, aes(x=Length_cm, y=Weight_kg)) +
  geom_point() +
  theme_sleek()
# weight_length

ggsave(filename="plots/assessment/weight_length.pdf", weight_length,
       width=150, height=100, units="mm", dpi=300)


# Maturity --------------------------------------------------------------------
# Age and maturity
age_maturity <- ggplot(maturity_table, aes(x=age, y=maturity)) +
  geom_point() +
  theme_sleek() 
# age_maturity

ggsave(filename="plots/assessment/age_maturity.pdf", age_maturity,
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
age_data <- read.csv("data/assessment/age_comp.csv")

catch_age_comp <- ggplot(age_data, aes(x=year, y=catch, fill=age)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_sleek() +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(~source)
# catch_age_comp

ggsave(filename="plots/assessment/catch_age_comp.pdf", catch_age_comp,
       width=300, height=100, units="mm", dpi=300)


# Catch rates -----------------------------------------------------------------
all_catch <- read.csv("data/assessment/catch_source.csv")

all_catch$source <- factor(all_catch$source, 
                           levels = c("freezer-trawler", "JV", "shoreside",
                                      "catcher-processor", "mothership", 
                                      "research", "unidentified"))

catch <- ggplot(all_catch, aes(y=catch, x=year, fill=source)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  theme_sleek() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~nation, ncol=2)

ggsave(filename="plots/assessment/catch.pdf", catch,
       width=300, height=150, units="mm", dpi=300)


# Survey ----------------------------------------------------------------------
survey_history <- read.csv(paste0(path, "survey-history.csv"))
survey_history$year <- as.factor(survey_history$year)

survey_biomass <- ggplot(survey_history, aes(y=biomass, x=year)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=biomass-cv, ymax=biomass+cv), width=.2) +
  theme_sleek()
# survey_biomass

ggsave(filename="plots/assessment/survey_biomass.pdf", survey_biomass,
       width=150, height=100, units="mm", dpi=300)

# Age - Length transition matrix ----------------------------------------------
# Pollock from Grant's CEATTLE implementation
test_data <- read.csv("data/assessment/age_trans.csv")

test2 <- melt(test_data, id.vars = "Age", 
              variable.name = "length_cat", value.name = "prop")

age_trans <- ggplot(test2, aes(x=length_cat, y=prop)) +
  geom_point() +
  theme_sleek() +
  facet_wrap(~Age, ncol=2)
age_trans
