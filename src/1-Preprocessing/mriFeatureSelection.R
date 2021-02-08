########################################################################
#
# mriFeatureSelection.R
#
# The following R script computes the MRI brain morphological features 
# more likely to separate healthy controls from Parkinson's Disease 
# (PD) patients. 
#
# To evaluate predictive power of MRI features we compute the adjusted 
# P-value of a linear model fitting MRI feature value and subject 
# status (the intercept is corrrected by age and gender of the 
# individual).
#
# The data used come from PPMI data portal and were extracted by 
# processing MRI images available on the database.
#
# TO DO: add very brief description of MRI
#
########################################################################


# R packages required to run the srcipt:
#   - tidyverse
#   - dplyr
#   - assertthat
if (!require("tidyverse", character.only = TRUE))
{
  install.packages("tidyverse", dependencies = TRUE)
  if (!require("tidyverse", character.only = TRUE))
  {
    stop("tidyverse package not found")
  }
}

if (!require("dplyr", character.only = TRUE))
{
  install.packages("dplyr", dependencies = TRUE)
  if (!require("dplyr", character.only = TRUE))
  {
    stop("dplyr package not found")
  }
}

if (!require("assertthat", character.only = TRUE))
{
  install.packages("assertthat", dependencies = TRUE)
  if (!require("assertthat", character.only = TRUE))
  {
    stop("assertthat package not found")
  }
}


# load the required packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(assertthat))


# function to compute the p-value of a linear model 
compute.lm <- function(feature)
{
  fit <- lm(
    as.formula(
      paste(
        feature, "~", paste("ENROLL_CAT", "eTIV", "Age", "GENDER", sep = "+"),
        sep = ""
      )
    ),
    data = eu.mri.baseline
  )
  
  p <- summary(fit)$coeff["ENROLL_CAT",4]
  
  return(p)
}


# retrieve data directory path
data.dir <- "~/Desktop/IGproject/IGPD/data/"  # modify accordingly

# get IDs of subjects with European ancestry retrieved during the QC step
eu.pats <- read.csv(
  paste(data.dir, "genotyping/PPMI_eu_woswedd_ds.fam", sep = ""),
  header = FALSE,
  sep = " "
)[,1:2]

# retrieve MRI imaging features
mri.feat <- read.csv(
  paste(data.dir, "imaging/MRI/mriFeatures.csv", sep = ""),
)[-1] # remove first column --> unuseful index

# get only data retrieved during the first visit (baseline)
mri.baseline <- mri.feat[mri.feat$EVENT_ID == "BL",] %>% distinct(PATNO, .keep_all = TRUE)

# get only data related to individuals with european ancestry
eu.mri.baseline <- mri.baseline[mri.baseline$PATNO %in% eu.pats$V1,]
if (!are_equal(nrow(eu.pats), nrow(eu.mri.baseline)))
    stop("wrong number of MRI features retrieved")


# load the enrolment cathegory information
pat.info <- read.csv(
  paste(data.dir, "patient_docs/Patient_Status.csv", sep = "")
) 
eu.pat.info <- pat.info[pat.info$PATNO %in% eu.pats$V1,]
if (!are_equal(nrow(eu.pats), nrow(eu.mri.baseline)))
  stop("wrong number of patients retrieved")
eu.mri.baseline$ENROLL_CAT <- ifelse(eu.pat.info$ENROLL_CAT == "PD", 1, 0)

# store the baseline MRI data for subjects with european ancestry
write.csv(
  eu.mri.baseline[,c(1,632:ncol(eu.mri.baseline))],
  file = paste(data.dir, "imaging/MRI/mriFeatures_eu_woswedd.csv", sep = ""),
  quote = FALSE,
  row.names = FALSE
)

# initialize result report table
features <- data.frame(
  Name = colnames(eu.mri.baseline)[632:(ncol(eu.mri.baseline) - 3)],
  p = rep(-1, length(colnames(eu.mri.baseline)[632:(ncol(eu.mri.baseline) - 3)])),
  p.adj = rep(-1, length(colnames(eu.mri.baseline)[632:(ncol(eu.mri.baseline) - 3)]))
  )

# compute lm p-values
p <- sapply(as.list(as.character(features$Name)), compute.lm)

features$p <- p
features <- features[!is.na(features$p),]
features$p.adj <- p.adjust(features$p, method = "BH")
features <- features[order(features$p.adj),]

# store the result in a file named mriFeatureRank.csv
write.csv(
  features,
  file = paste(data.dir, "imaging/MRI/mriFeaturesRank.csv", sep = ""),
  quote = FALSE,
  row.names = FALSE
)

