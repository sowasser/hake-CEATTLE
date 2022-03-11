# Script for any manipulation of arrowtooth flounder assessment data needed for
# inclusion in CEATTLE

library(dplyr)
library(ggplot2)

# Consolodate selectivity data
sel_initial <- read.csv("data/assessment/Arrowtooth_Flounder/selectivity.csv")

# sel <- sel_initial %>% group_by(Fleet, Year, Sex) %>%
#   summarize()

# Calculate sex ratio ---------------------------------------------------------
ages <- 2:35
ratio <- exp((0.166-0.274) * ages)
print(ratio)

# Plot to check ratios against assessment doc
sex_ratios <- as.data.frame(cbind(ages, ratio))
ggplot(sex_ratios, aes(x=ages, y=ratio)) +
  geom_line() +
  theme_sleek()

write.csv(sex_ratios, "data/assessment/Arrowtooth_Flounder/ATF_sexratio.csv")
