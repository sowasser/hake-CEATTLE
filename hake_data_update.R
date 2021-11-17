# Script for data manipulation for adapting the hake stock assessment for 
# input into CEATTLE

library(reshape2)
library(dplyr)

path <- "~/Desktop/Local/hake-assessment-master/data/"

# Catch data ------------------------------------------------------------------
# Comes from multiple sources in the US and Canada

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
us_unid <- as.data.frame(cbind(year = us_unid[, 3], 
                               nation = rep("US", length(us_unid[, 3])), 
                               source = rep("unidentified", length(us_unid[, 3])),
                               month = us_unid[, 2], 
                               catch = us_unid[, 4]))

monthly_catch <- rbind(catch_data("can-ft-catch-by-month.csv", "Canada", "freezer-trawler"),
                       catch_data("can-jv-catch-by-month.csv", "Canada", "JV"),
                       catch_data("can-ss-catch-by-month.csv", "Canada", "shoreside"),
                       catch_data("us-cp-catch-by-month.csv", "US", "catcher-processor"),
                       catch_data("us-ms-catch-by-month.csv", "US", "mothership"),
                       catch_data("us-research-catch-by-month.csv", "US", "research"),
                       catch_data("us-shore-catch-by-month.csv", "US", "shoreside"),
                       us_unid)

monthly_catch$catch <- as.numeric(monthly_catch$catch)
all_catch <- monthly_catch %>% group_by(year, nation, source) %>%
  summarize(total_catch = sum(catch))

write.csv(all_catch, "data/all_catch.csv")
