# Exploring the hake dynamics from the stock assessment data:
# https://github.com/pacific-hake/hake-assessment/tree/master/data

library(ggplot2)
library(ggsidekick)

maturity <- read.csv("~/Desktop/Local/hake-assessment-master/data/maturity-table.csv")

ggplot(maturity, aes(x = age, y = avg.wt)) +
  geom_point() +
  theme_sleek()