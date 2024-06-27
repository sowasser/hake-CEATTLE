#' Script for adapting the Stock Synthesis 3 files (and tables from the document) 
#' for the 2024 hake assessment found here 
#' (https://github.com/pacific-hake/m2/tree/main/inst/extdata/ss3) and here 
#' (https://github.com/pacific-hake/hake-assessment/tree/master/data-tables)
#' to the format needed for the Rceattle input.

library(here)
library(dplyr)
library(reshape2)

# Create spawn weight (maturity * weight-at-age) ------------------------------
# Read in weight-at-age table pulled from SS3 input
wtatage <- read.table(here("data", "assessment", "2024", "wtatage.txt"))

# Make labels for columns
a <- 1:20
for(i in 1:length(a)) {
  a[i] <- paste0("a", i)
}
a <- append("a0", a)
colnames(wtatage) <- c("Yr", "Seas", "Sex", "Bio_Pattern", "BirthSeas", "Fleet", a)

# Select fishery weight-at-age for the years needed & get rid of extra columns
fishery_wt <- wtatage %>% filter(Fleet == 1 & Yr %in% 1975:2023) 
fishery_wt <- fishery_wt[, c(1, 7:27)]

# Read in maturity ogives from assessment document & filter for "null" model (at least for now).
ogives <- read.csv(here("data", "assessment", "2024", "maturity-ogives.csv")) %>%
  filter(model == "Null") %>%
  select(!model)

# Reshape to wide format to match weight-at-age
maturity <- dcast(ogives, formula = year ~ age)
colnames(maturity) <- c("Yr", a[1:16])

# Add values for years previous to 2009 & remove projection years after 2023
for(i in 1975:2008) {
  x <- c(i, maturity[1, 2:17])
  maturity[nrow(maturity) +1, ] <- x
}
maturity <- maturity[order(maturity$Yr), ] %>%
  filter(Yr %in% 1975:2023)

# Add columns for a16-20 (with a15 values)
maturity$a16 <- maturity[, 17]
maturity$a17 <- maturity[, 17]
maturity$a18 <- maturity[, 17]
maturity$a19 <- maturity[, 17]
maturity$a20 <- maturity[, 17]

# Multiply together to get spawning weight & write .csv
spawn_wt <- fishery_wt[, -1] * maturity[, -1]

spawn_wt <- cbind.data.frame(Year = maturity$Yr, spawn_wt)
write.csv(spawn_wt, file = here("data", "assessment", "2024", "spawn_wtatage.csv"))


# Age transition matrix / age-length key --------------------------------------
# Read in larger maturity table w/ length-at-age from assessment document
matdat <- read.csv(here("data", "assessment", "2024", "maturity-samples.csv"))

# Set up bins 
hake_lbin <- c(0, seq(20, 70, by = 2), 999)
hake_age_bin <- c(0:14 + 0.5, 99)
hake_ages <- 1:20

# Cut into bins
matdat$BIN <- cut(matdat$length_cm, breaks = hake_lbin)
levels(matdat$BIN) <- 1:length(hake_lbin[-1])

matdat$AgeBIN <- cut(matdat$age, breaks = hake_age_bin)
levels(matdat$AgeBIN) <- hake_ages

# Set up age-length key
mat_new <- matdat[-which(is.na(matdat$AgeBIN)),]
mat_table <- with(mat_new,table(BIN,AgeBIN))
alk <- prop.table(mat_table, margin=1)
alk[is.na(alk)] <- 0

# Save as age_trans_matrix in CEATTLE input excel sheet
age_trans_matrix <- t(as.data.frame.matrix(alk))
write.csv(age_trans_matrix, here("data", "assessment", "2024", "age_trans_matrix.csv"), row.names = FALSE)
