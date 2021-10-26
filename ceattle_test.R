# Install necessary packages --------------------------------------------------
# install.packages("devtools")

# Install TMB and rtools https://cran.r-project.org/bin/windows/Rtools/
# Instructions here for non-pc: https://github.com/kaskr/adcomp/wiki/Download
# install.packages('TMB', type = 'source')
TMB::runExample(all = TRUE)  # see if TMB works

# Install Rceattle - https://github.com/grantdadams/Rceattle
# devtools::install_github("grantdadams/Rceattle")
library(Rceattle)

