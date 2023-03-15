### Age-length key method for ageing ------------------------------------------
# Set up hake age-length key
hake_lbin <- c(0, seq(20, 70, by = 2), 999)
hake_age_bin <- c(0:14 + 0.5, 99)
hake_ages <- 1:15

maturity <- read.csv("~/Desktop/Local/hake-CEATTLE/Resources/hake-assessment-master/data/hake-maturity-data.csv")
maturity$BIN <- cut(maturity$Length_cm, breaks = hake_lbin)
levels(maturity$BIN) <- 1:length(hake_lbin[-1])

maturity$AgeBIN <- cut(maturity$Age, breaks = hake_age_bin)
levels(maturity$AgeBIN) <- hake_ages

hake <- maturity[-which(is.na(maturity$AgeBIN)),]
hake_table <- with(hake,table(BIN,AgeBIN))
alk_hake <- prop.table(hake_table, margin=1)
alk_hake[is.na(alk_hake)] <- 0

# Save for age_trans_matrix in CEATTLE input excel sheet
age_trans_matrix <- t(as.data.frame.matrix(alk_hake))
write.csv(age_trans_matrix, "data/assessment/age_trans_matrix.csv", row.names = FALSE)
