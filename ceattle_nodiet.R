# Run CEATTLE for hake data with no diet data 

library(Rceattle)
library(readxl)

# Read in data from excel file ------------------------------------------------
read_excel_data <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

hake_nodiet <- read_excel_data("data/input_nodiet.xlsx")