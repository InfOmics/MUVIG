########################################################################
#
# plot_datscan_gwas.R
#
# The following R script plots the results obtained during Genome wide
# Association Study using DaTSCAN uptake values to improve SNP-phenotype
# associations.
#
# We will compute, for each DaTSCAN feature considered (CAUDATE_R, 
# CAUDATE_L, PUTAMEN_R, PUTAMEN_L) the corresponding manhattan plot,
# showing the statistical significance of SNP-phenotype association for 
# each genetic variants.
#
# We will also compute the quantile-quantile plots (qq plots) for each 
# GWAS analysis results.
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
source("Manhattan_plot.R") # contains the function to generate the Manhattan plot

genotyping.dir <- "../../data/genotyping/Individual_View_Datscan/"
results.dir    <- "../../results/individualView/DaTSCAN/"
phenotype.prefix <- paste(genotyping.dir, "indview_datscan.", sep = "")

# phenotypes
phenotypes = c("CAUDATE_R","CAUDATE_L","PUTAMEN_R","PUTAMEN_L")

# compute manhanttan and qq plots for each DaTSCAN feature
for (pheno in phenotypes){
  snp.pheno.ass = read.table(
    paste(paste(phenotype.prefix, pheno, sep = ""), ".assoc.linear", sep = ""),
    header = TRUE
  )
  
  out.fn <- paste(paste(results.dir, pheno, sep = ""), "_Manhattan.png", sep = "")
  out.fn.qq <- paste(paste(results.dir, pheno, sep = ""), "_QQ.png", sep = "")

  # manhanttan plot
  ManhattanGenerator(snp.pheno.ass,out.fn,pheno)
  
  # qq plot
  png(out.fn.qq, width = 24, height = 16, units = "in", res = 275) 
  qq(snp.pheno.ass$P, main = sprintf("%s - QQ plot", pheno), col = "blue4")
  dev.off()  # close ostream
}









