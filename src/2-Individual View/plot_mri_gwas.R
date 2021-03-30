########################################################################
#
# plot_mri_gwas.R
#
# The following R script plots the results obtained during Genome wide
# Association Study using MRI structural measures to improve SNP-phenotype
# associations.
#
# We will compute, for each MRI feature considered 
# (rh_parahippocampal_volume, lh_parahippocampal_volume, 
# rh_parahippocampal_area, lh_parahippocampal_area) the corresponding
# manhattan plot, which shows the statistical significance of SNP-phenotype 
# association for each genetic variants.
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


genotyping.dir <- "../../data/genotyping/"
results.dir <-"../../results/individualView/MRI/"
phenotype.prefix <- paste(genotyping.dir, "Individual_View_MRI/indview_mri.", sep = "")

# phenotypes
phenotypes = c(
  "rh_parahippocampal_volume",
  "lh_parahippocampal_volume",
  "rh_parahippocampal_area",
  "lh_parahippocampal_area",
  "rh_parahippocampal_thickness",
  "lh_parahippocampal_thickness"
)

# compute manhanttan and qq plots for each DaTSCAN feature
for (pheno in phenotypes){
  snp.pheno.ass = read.table(
    paste(paste(phenotype.prefix, pheno, sep = ""), ".assoc.linear", sep = ""),
    header = TRUE
  )
  
  out.fn <- paste(paste(results.dir, pheno, sep = ""), "_Manhattan.png", sep = "")
  out.fn.qq <- paste(paste(results.dir, pheno, sep = ""), "_QQ.png", sep = "")
  
  # manhanttan plot
  ManhattanGenerator(snp.pheno.ass, out.fn, pheno)
  
  # qq plot
  png(out.fn.qq, width = 24, height = 16, units = "in", res = 275) 
  qq(snp.pheno.ass$P, main = sprintf("%s - QQ plot", pheno), col = "blue4")
  dev.off()  # close ostream
}
