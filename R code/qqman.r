setwd("/Users/guglielmo/Desktop/Università/Magistrale/Tesi/RECOMPUTE/GWAS2/MRI")
library(qqman)
check=function(x) tryCatch(if(class(x) == 'logical') 1 else 1, error=function(e) 0) 


name = c("CAUDATE_L","CAUDATE_R","PUTAMEN_R","PUTAMEN_L")
name2 = c("rh_parahippocampal_area","lh_parahippocampal_area",
          "rh_parahippocampal_volume","lh_parahippocampal_volume")

  
for (i in name2){
  if(check(name)){
    gwas = read.csv(sprintf("GWAS_DATScan_%s.csv", i))
  } else{
    gwas = read.csv(sprintf("GWAS_MRI_%s_Without_SWEDD_win.csv", i))  
  }

png(sprintf("PLOT/%s_Without_Swedd_win.png", i), width = 24, height = 16, units = "in", res = 275)
par(mfrow=c(2,1)) 
manhattan(gwas, main = sprintf("Manhattan Plot of %s", i), suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), 
          chr="CHR", bp="BP", p="P",ylim = c(0, 8), col = c("red2", "green3") , annotatePval = 0.01)
qq(gwas$P, main = sprintf("Q-Q plot of %s", i), col = "blue4")
dev.off()
}
