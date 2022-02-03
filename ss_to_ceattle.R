# Script for adapting Stock Synthesis (SS3) outputs for CEATTLE

wt <- read.table("data/assessment/wtatage.txt")
names <- toString(0:20)
colnames(wt) <- c("yr", "seas", "sex", "bio_pattern", "birthseas", "fleet", 
                  "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", 
                  "a10", "a11", "a12", "a13", "a14", "a15", "a16", "a17", 
                  "a18", "a19", "a20")
wt2 <- wt[order(wt$fleet), ]
