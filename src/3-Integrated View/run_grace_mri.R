########################################################################
#
# run_grace_mri.R
#
# The following R script computes growth curve models fitting DATScan 
# measure evolution among the European subjects of PPMI.
#
# The models are computed using GRACE tool. The growth curve are built
# to fit the DATScan measure changes as functions of subjects age.
# 
# After GRACE computation we have a single measure describing the 
# differences in DATScan uptake values among the PPMI subjects.
#
# This is the fundamental step of the Integrated View phase, since 
# during this step we merge different (but closely related) phenotypic
# values to a single unifying measure (gamma0).
#
# The new values are used as single phenotypic information to use 
# during GWAS computation
#
########################################################################


# R packages required to run the script:
#   - devtools
#   - grace
#   - fda
#   - mvtnorm
#   - tidyverse
if(!require("devtools", character.only = TRUE))
{
  install.packages("devtools")
  if(!require("devtools", character.only = TRUE))
  {
    stop("devtools package not found")
  }
}

if(!require("grace", character.only = TRUE))
{
  library(devtools)
  install_bitbucket("grace", "mdonohue")
  if(!require("grace", character.only = TRUE))
  {
    stop("grace package not found")
  }
}

if(!require("fda", character.only = TRUE))
{
  install.packages("fda")
  if(!require("fda", character.only = TRUE))
  {
    stop("fda package not found")
  }
}

if(!require("mvtnorm", character.only = ))
{
  install.packages("mvtnorm")
  if(!require("mvtnorm", character.only = TRUE))
  {
    stop("mvtnorm package not found")
  }
}

if(!require("tidyverse", character.only = TRUE))
{ 
  install.packages("tidyverse")
  if(!require("tidyverse", character.only = TRUE))
  {
    stop("tidyverse package not found")
  }
}

# load the required packages
library(devtools)
library(grace) 
library(fda) 
library(mvtnorm) 
library(tidyverse)
library(ggplot2) 


options(stringsAsFactors = FALSE) 

# load DATScan measures
mri.dir <- "../../data/imaging/MRI/"
mri <- read.csv(paste(mri.dir, "MRIFeatures_eu_woswedd.csv", sep = ""))
mri <- mri[, c(1,3:267)]
mri <- mri[complete.cases(mri),]

# load phenotype information (computed during prev step)
pheno.fn <- "../../data/genotyping/phenotype_mri.txt"
pheno <- read.csv(pheno.fn, sep = " ")

# load covariate information (computed during prev step)
covar.fn = "../../data/genotyping/covariate_mri.txt"
covar <- read.csv(covar.fn, sep = "\t")

# create df for grace run
features <- colnames(mri)[2:length(colnames(mri))]
# for (f in features)
# {
#   mri[f] <- (mri[[f]] - min(mri[[f]], na.rm = T)) / sd(mri[[f]], na.rm = T) 
# }
grace.df <- data.frame(
  id = rep(mri$PATNO, each = 265),  # 265 features
  argvals = rep(covar[covar$IID %in% mri$PATNO,]$AGE, each = 265),
  group = as.factor(rep(pheno[pheno$IID %in% mri$PATNO,]$ENROLL_CAT, each = 265)),
  Y = as.numeric(as.list(t(mri[, 2:ncol(mri)]))),
  outcome = rep(features, length(mri$PATNO))
)

# run grace
grace.fits <- grace(id = grace.df$id,
                    argvals = grace.df$argvals,
                    y = grace.df$Y,
                    outcome = grace.df$outcome,
                    group = grace.df$group,
                    plots = T)

# store grace results
for (f in features)
{
  out.name <- paste(datscan.dir, f, "_grace.csv", sep = " ")
  write.csv(x = grace.fits$fits[[f]]$subset, file = out.name, quote = FALSE, row.names = FALSE) 
}

