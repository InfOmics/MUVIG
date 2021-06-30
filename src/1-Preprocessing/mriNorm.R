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

mri.dir <- "../../data/patient_data/"

# read the DaTSCAN 
mri.ceu <- read.csv(
  paste(mri.dir, "MRI_CEU_filt.csv", sep = "")
)

# normalize values with rank based inverse normal transform
mri.ceu.norm <- mri.ceu
mri.ceu.norm[,2:ncol(mri.ceu.norm)] <- apply(
  mri.ceu.norm[,2:ncol(mri.ceu.norm)], 2, RankNorm
)

# write the results
write.csv(
  mri.ceu.norm,
  file = paste(mri.dir, "MRI_CEU_norm.csv", sep = ""),
  quote = FALSE,
  row.names = FALSE
)





