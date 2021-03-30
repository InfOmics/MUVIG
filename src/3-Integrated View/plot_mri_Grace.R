########################################################################
#
# plot_DATScan_GRACE.R
#
# The following R script plots the results obtained after GRACE results
#
# We will compute the manhattan plot and the quantile-quantile (qq) plot.
#
########################################################################


# R packages required to run the script:
#   - qqman
if(!require("qqman", character.only = TRUE))
{
  install.packages("qqman")
  if(!require("qqman", character.only = TRUE))
  {
    stop("qqman package not found")
  }
}


# load required packages
suppressPackageStartupMessages(library(qqman))
source("../2-Individial View/Manhattan_plot.R")

data.dir    = "../../data/genotyping/Integrated_View_MRI"
results.dir = "../../results/integratedView/MRI"

gwas = read.table(paste0(data.dir,"/Grace_MRI.assoc.linear"), header = T)

out.fn    = paste0(results.dir,"/Manhattan_mri_Grace.png")
out.fn.qq = paste0(results.dir,"/QQ_mri_Grace.png")

ManhattanGenerator(gwas, out.fn, "GRACE")

png(out.fn.qq, width = 24, height = 16, units = "in", res = 275)
qq(gwas$P, main = "Q-Q plot", col = "blue4")
dev.off()

