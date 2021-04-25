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
# processing MRI images available on the data portal.
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
    data = mri.foi.bl
  )
  
  p <- summary(fit)$coeff["ENROLL_CAT",4]
  
  return(p)
}


# retrieve data directory path
data.dir <- "~/Desktop/MUVIG-main/data/"  # modify accordingly

# retrieve MRI imaging features
mri.foi.bl <- read.csv(
  paste(data.dir, "imaging/MRI/mriFeatures_foi_bl.csv", sep = ""),
)

# load the enrolment cathegory information
pat.info <- read.csv(
  paste(data.dir, "patient_docs/Patient_Status.csv", sep = "")
) 
eu.pat.info <- pat.info[pat.info$PATNO %in% mri.foi.bl$PATNO,]
if (!are_equal(nrow(eu.pat.info), nrow(mri.foi.bl)))
  stop("wrong number of patients retrieved")
mri.foi.bl$ENROLL_CAT <- ifelse(eu.pat.info$ENROLL_CAT == "PD", 2, 1)

# initialize result report table
features <- data.frame(
  Name = colnames(mri.foi.bl)[6:(ncol(mri.foi.bl) - 1)],
  p = rep(-1, length(colnames(mri.foi.bl)[6:(ncol(mri.foi.bl) - 1)])),
  p.adj = rep(colnames(mri.foi.bl)[6:(ncol(mri.foi.bl) - 1)])
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

