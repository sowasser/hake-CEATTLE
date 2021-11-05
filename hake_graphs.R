# Exploring the hake dynamics from the stock assessment data:
# https://github.com/pacific-hake/hake-assessment/tree/master/data

library(ggplot2)
library(ggsidekick)


# Weight & length at age ------------------------------------------------------
maturity_table <- read.csv("~/Desktop/Local/hake-assessment-master/data/maturity-table.csv")
hake_maturity_data <- read.csv("~/Desktop/Local/hake-assessment-master/data/hake-maturity-data.csv")

weight_age <- ggplot(maturity_table, aes(x=age, y=avg.wt)) +
  geom_point() +
  theme_sleek()
# weight_age

ggsave(filename="plots/weight_age.pdf", weight_age,
       width=150, height=100, units="mm", dpi=300)

weight_age_yearly <- ggplot(hake_maturity_data, aes(x=Age, y=Weight_kg)) +
  geom_point() +
  theme_sleek() +
  facet_wrap(~Year, ncol = 3) 
# weight_age_yearly

ggsave(filename="plots/weight_age_yearly.pdf", weight_age_yearly,
       width=300, height=200, units="mm", dpi=300)

length_age_yearly <- ggplot(hake_maturity_data, aes(x=Age, y=Length_cm)) +
  geom_point() +
  theme_sleek() +
  facet_wrap(~Year, ncol = 3) 
# length_age_yearly

ggsave(filename="plots/length_age_yearly.pdf", length_age_yearly,
       width=300, height=200, units="mm", dpi=300)
