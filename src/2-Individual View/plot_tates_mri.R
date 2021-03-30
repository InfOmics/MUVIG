########################################################################
#
# plot_mri_tates.R
#
# The following R script plots the results obtained after single-trait
# based statistical significance of GWAS results made with TATES.
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
source("Manhattan_plot.R")

genotyping.dir <- "../../data/genotyping/"
tates.dir <- paste(genotyping.dir, "tates_mri_wd/", sep = "")
results.dir <-"../../results/individualView/MRI/"

# compute manhanttan and qq plot
snp.pheno.st.ass = read.table(
  paste(tates.dir, "tates_mri_results", sep = ""),
)
colnames(snp.pheno.st.ass) <- c("CHR", "SNP", "PHENO", "P")
snp.pheno.st.ass = cbind(
  snp.pheno.st.ass, data.frame(
    BP = read.table(
      paste(
        paste(
          paste(genotyping.dir, "Individual_View_MRI/indview_mri.", sep = ""), 
          "rh_parahippocampal_volume", sep = ""
        ), 
        ".assoc.linear", sep = ""
      ), header = T
    )$BP
  )
)

out.fn    <- paste(paste(results.dir, "tates_st_correction_Manhattan", sep = ""), ".png", sep = "")
out.fn.qq <- paste(paste(results.dir, "tates_st_correction_QQ", sep = ""), ".png", sep = "")

# manhanttan plot
ManhattanGenerator(snp.pheno.st.ass, 
                   out.fn, 
                   "TATES ST correction")

# qq plot
png(out.fn.qq, width = 24, height = 16, units = "in", res = 275)
qq(snp.pheno.st.ass$P, main = "TATES ST correction - QQ plot", col = "blue4")
dev.off()  # close ostream  
