# Run CEATTLE for hake data with no diet data 

library(Rceattle)

mydata <- Rceattle::read_data( file = "data/input_nodiet.xlsx")
mydata$est_M1 <- c(0,0,0)
