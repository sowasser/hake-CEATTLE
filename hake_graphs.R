# Exploring the hake dynamics from the stock assessment data:
# https://github.com/pacific-hake/hake-assessment/tree/master/data

library(ggplot2)
library(ggsidekick)
library(reshape2)
library(viridis)


# Weight & length at age ------------------------------------------------------
maturity_table <- read.csv("~/Desktop/Local/hake-assessment-master/data/maturity-table.csv")
hake_maturity_data <- read.csv("~/Desktop/Local/hake-assessment-master/data/hake-maturity-data.csv")

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


# Age composition -------------------------------------------------------------
us_cp_age_data <- read.csv("~/Desktop/Local/hake-assessment-master/data/us-cp-age-data.csv")
us_ms_age_data <- read.csv("~/Desktop/Local/hake-assessment-master/data/us-ms-age-data.csv")
us_shore_age_data <- read.csv("~/Desktop/Local/hake-assessment-master/data/us-shore-age-data.csv")

us_cp_age_data <- cbind(us_cp_age_data[, -3], source=rep("cp", length(us_cp_age_data[, 1])))
us_ms_age_data <- cbind(us_ms_age_data[, -3], source=rep("ms", length(us_ms_age_data[, 1])))
us_shore_age_data <- cbind(us_shore_age_data[, -3], source=rep("shore", length(us_shore_age_data[, 1])))
age_data_wide <- rbind(us_cp_age_data, us_ms_age_data, us_shore_age_data)

age_data <- melt(age_data_wide, id.vars = c("year", "n.fish", "source"))
colnames(age_data) <- c("year", "n.fish", "source", "age", "proportion")
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
us_catch_by_month <- read.csv("~/Desktop/Local/hake-assessment-master/data/us-catch-by-month.csv")

monthly_catch <- ggplot(us_catch_by_month, aes(y=Catch.MT, x=month)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_sleek() +
  facet_wrap(~year, ncol=3)
# monthly_catch

ggsave(filename="plots/monthly_catch.pdf", monthly_catch,
       width=300, height=300, units="mm", dpi=300)
