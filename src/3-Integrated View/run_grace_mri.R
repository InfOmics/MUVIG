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
suppressPackageStartupMessages( c( library(devtools),
                                   library(grace),
                                   library(fda),
                                   library(mvtnorm),
                                   library(tidyverse) ))

options(stringsAsFactors = FALSE) 

genotyping.dir = "../../data/genotyping/"
indview.mri    = paste0(genotyping.dir, "Individual_View_MRI/")
grace.dir      = "Integrated_View_MRI"
results.dir    = "../../results/"

# create the folder that will contains the results
dir.create(file.path(genotyping.dir, grace.dir), showWarnings = TRUE)

# load phenotype information (computed during prev step)
pheno = read.csv(paste0(indview.mri,"phenotype_mri.txt"), sep = " ")

# to run properly grace we need to remove all NA values in the phenotype file
pheno = pheno[complete.cases(pheno),]

# load covariate information (computed during prev step)
covar = read.csv(paste0(indview.mri,"covariate_mri.txt"), sep = " ")

# we need to remove the patients also in the covariate file to have the same patients set
covar = covar[covar$FID %in% pheno$FID,]

# create df for grace run
features <- colnames(pheno)[4:length(colnames(pheno))]

times = length(features)

grace.df <- data.frame(
  id = rep(pheno$FID, each = times),      # 6 features
  argvals = rep(covar$AGE, each = times),
  group = as.factor(rep(pheno$ENROLL_CAT, each = times)),
  Y = as.numeric(as.list(t(pheno[, features]))),
  outcome = rep(features, length(pheno$FID))
)

# before run grace is suggested to move to another folder in order to store the resulting plot
dir.create(file.path(results.dir, "integratedView"), showWarnings = TRUE)
dir.create(file.path(results.dir, "integratedView/MRI"), showWarnings = TRUE)

setwd(paste0(results.dir, "integratedView/MRI"))

grace.fits <- grace(id = grace.df$id,
                    argvals = grace.df$argvals,
                    y = grace.df$Y,
                    outcome = grace.df$outcome,
                    group = grace.df$group,
                    plots = T)

setwd("../")

# store grace table results
for (f in features){
  out.name <- paste0(f, "_grace.csv")
  path = paste0(genotyping.dir, grace.dir, "/", out.name)
  write.csv(x = grace.fits$fits[[f]]$subset, file = path, quote = FALSE, row.names = FALSE) 
}
