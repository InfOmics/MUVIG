########################################################################
#
# mriNorm.R
#
# The following R script normalizes the MRI structural values using 
# rank-based inverse normal transform (r-INT) procedure.
#
# We chose r-INT because it is less susceptible to outlier values during
# normalization.
#
# To avoid possible errors and biases due to outliers and different 
# value ranges among MRI features when associating SNP with phenotypes, 
# we need to normalize the values.
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

mri.dir <- "../../data/patients_data/"

# read the DaTSCAN 
mri.parahippo <- read.csv(
  paste(mri.dir, "MRI_parahippo.csv", sep = "")
)

# normalize values with rank based inverse normal transform
mri.parahippo.norm <- mri.parahippo
mri.parahippo.norm[,5:(ncol(mri.parahippo.norm) - 1)] <- apply(
  mri.parahippo.norm[,5:(ncol(mri.parahippo.norm) - 1)], 2, RankNorm
)

# write the results
write.csv(
  mri.parahippo.norm,
  file = paste(mri.dir, "MRI_parahippo_norm.csv", sep = ""),
  quote = FALSE,
  row.names = FALSE
)





