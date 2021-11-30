# Constructing what CEATTLE calls an age transition matrix and whay fishR calls
# Age-Length Keys - essentially the proportional lengths of fish at each age.
# Vignette here: http://derekogle.com/fishR/examples/oldFishRVignettes/AgeLengthKey.pdf

library(dplyr)
library(FSA)

path <- "~/Desktop/Local/hake-assessment-master/data/"

# Extract length at age from assessment data 
hake_maturity_data <- read.csv(paste0(path, "hake-maturity-data.csv"))
length_age_all <- na.omit(hake_maturity_data[, 11:12])

# Bin all ages 15+ into age 15 (following other parts of the assessment)
length_age_all$Age[length_age_all$Age >= 15] <- 15

# Check number of unique lengths for binning
# length(unique(length_age_all$Length_cm))

# Determine length categories - set w to a value that will give # lengths = ages
len.age.cat <- lencat(~Length_cm, data=length_age_all, startcat=16, w=4.2)

# Create length-at-age matrix, then find proportions
len.age.raw <- with(len.age.cat, table(LCat, Age))
len.age.prop <- prop.table(len.age.raw, margin = 1)

# Transpose to match CEATTLE and export
len.age <- as.data.frame.matrix(t(len.age.prop))
colnames(len.age) <- c("Length1", "Length2", "Length3", "Length4", "Length5", 
                       "Length6", "Length7", "Length8", "Length9", "Length10",
                       "Length11", "Length12", "Length13", "Length14", "Length15")

write.csv(len.age, "data/length_age_matrix.csv", row.names = FALSE)

