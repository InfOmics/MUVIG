########################################################################
#
# datscanNorm.R
#
# The following R script normalizes the DaTSCAN uptake values using 
# rank-based inverse normal transform (r-INT) procedure.
#
# We chose r-INT because it is less susceptible to outlier values during
# normalization.
# 
# In our study we had 4 DaTSCAN features:
#   - caudate right uptake values 
#   - caudate left uptake values
#   - putamen right uptake values
#   - putamen left uptake values
#
# The uptake measurements are continous values, but they have different
# ranges of values (putamen != caudate range).
#
# To avoid possible errors and biases due to the different value ranges 
# when associating SNP with phenotypes (PD in our study), we need to
# normalize the values.
#
########################################################################


# R packages required to run the script:
#   - RNOmni
if(!require("RNOmni", character.only = TRUE))
{
  install.packages("RNOmni")
  if(!require("RNOmni", character.only = TRUE))
  {
    stop("RNOmni package not found")
  }
}


# load the required packages
suppressPackageStartupMessages(library(RNOmni))

datscan.dir <- "../../data/imaging/DaTSCAN/"

# read the DaTSCAN 
datscan.eu <- read.csv(
  paste(datscan.dir, "DATScan_Analysis_eu_fv.csv", sep = "")
)

# normalize values with rank based inverse normal transform
datscan.eu.norm <- datscan.eu
datscan.eu.norm[,4:ncol(datscan.eu.norm)] <- apply(
  datscan.eu.norm[,4:ncol(datscan.eu.norm)], 2, RankNorm
)

# write the results
write.csv(
  datscan.eu.norm,
  file = paste(datscan.dir, "DATScan_Analysis_eu_fv_norm.csv", sep = ""),
  quote = FALSE,
  row.names = FALSE
)







