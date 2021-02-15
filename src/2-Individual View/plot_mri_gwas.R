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

genotyping.dir <- "../../data/genotyping/"
results.dir <-"../../results/individualView/MRI/"
phenotype.prefix <- paste(genotyping.dir, "indview_mri_", sep = "")

# phenotypes
phenotypes = c(
  "rh_parahippocampal_volume",
  "lh_parahippocampal_volume",
  "rh_parahippocampal_area",
  "lh_parahippocampal_area"
)

# compute manhanttan and qq plots for each DaTSCAN feature
for (pheno in phenotypes){
  snp.pheno.ass = read.csv(
    paste(paste(phenotype.prefix, pheno, sep = ""), ".csv", sep = "")
  )
  
  out.fn <- paste(paste(results.dir, pheno, sep = ""), ".png", sep = "")
  png(out.fn, width = 24, height = 16, units = "in", res = 275) 
  
  par(mfrow=c(2,1)) 
  
  # manhanttan plot
  manhattan(
    snp.pheno.ass, main = sprintf("%s - Manhanttan plot", pheno), chr="CHR", bp="BP", p="P",
    ylim = c(0, 8), col = c("blue4", "orange3"), annotatePval = 0.005
  )
  
  # qq plot
  qq(snp.pheno.ass$P, main = sprintf("%s - QQ plot", pheno), col = "blue4")
  
  dev.off()  # close ostream
}
