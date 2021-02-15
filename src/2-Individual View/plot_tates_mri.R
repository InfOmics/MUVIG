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
    BP = read.csv(
      paste(
        paste(
          paste(genotyping.dir, "indview_mri_", sep = ""), 
          "rh_parahippocampal_volume", sep = ""
        ), 
        ".csv", sep = ""
      )
    )$BP
  )
)

out.fn <- paste(paste(results.dir, "tates_st_correction", sep = ""), ".png", sep = "")
png(out.fn, width = 24, height = 16, units = "in", res = 275) 

par(mfrow=c(2,1)) 

# manhanttan plot
manhattan(
  snp.pheno.st.ass, main = "TATES ST correction - Manhanttan plot", chr="CHR", bp="BP", p="P",
  ylim = c(0, 8), col = c("blue4", "orange3"), annotatePval = 0.005
)

# qq plot
qq(snp.pheno.st.ass$P, main = "TATES ST correction - QQ plot", col = "blue4")

dev.off()  # close ostream

