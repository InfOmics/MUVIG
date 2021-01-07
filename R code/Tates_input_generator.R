setwd("/Users/guglielmo/Desktop/Università/Magistrale/Tesi/RECOMPUTE/TATES")

answer <- readline(prompt="Do you want to use SWEDD patients: ")
MRI = 1

# ************************ FOR MRI ********************
if(MRI==1){
if (answer=="yes" || answer=="Yes") {
# with swedd
para_area_l = read.csv("../GWAS2/MRI/GWAS_MRI_lh_parahippocampal_area_With_SWEDD_win.csv")
para_area_r = read.csv("../GWAS2/MRI/GWAS_MRI_rh_parahippocampal_area_With_SWEDD_win.csv")
para_volume_r = read.csv("../GWAS2/MRI/GWAS_MRI_rh_parahippocampal_volume_With_SWEDD_win.csv")
para_volume_l = read.csv("../GWAS2/MRI/GWAS_MRI_lh_parahippocampal_volume_With_SWEDD_win.csv")
name = "With SWEDD"
} else {
# without swedd
para_area_l = read.csv("../GWAS2/MRI/GWAS_MRI_lh_parahippocampal_area_Without_SWEDD_win.csv")
para_area_r = read.csv("../GWAS2/MRI/GWAS_MRI_rh_parahippocampal_area_Without_SWEDD_win.csv")
para_volume_r = read.csv("../GWAS2/MRI/GWAS_MRI_rh_parahippocampal_volume_Without_SWEDD_win.csv")
para_volume_l = read.csv("../GWAS2/MRI/GWAS_MRI_lh_parahippocampal_volume_Without_SWEDD_win.csv")
name = "Without SWEDD"
  }
} else {
    # without swedd
    caudate_r   = read.csv("../GWAS2/GWAS_DATScan_CAUDATE_R.csv")
    caudate_l   = read.csv("../GWAS2/GWAS_DATScan_CAUDATE_L.csv")
    putamen_r   = read.csv("../GWAS2/GWAS_DATScan_PUTAMEN_R.csv")
    putamen_l   = read.csv("../GWAS2/GWAS_DATScan_PUTAMEN_L.csv")
}

#Create the File of p-values
if(MRI==1){
pvals = data.frame(para_area_l[c(1,2,9)],para_area_r$P,para_volume_r$P, para_volume_l$P)
colnames(pvals)[3:6] = c("P1-PAR","P2-PAR","P3-PVR","P4-PVL")
write.table(pvals,paste("MRI Results/WINSORIZE/",name,"/Pvals",sep=""), sep=" ",
            quote = F, row.names = F, col.names = F)
} else {
  name = ""
  pvals = data.frame(caudate_r[c(1,2,9)],caudate_l$P,putamen_r$P, putamen_l$P)
  colnames(pvals)[3:6] = c("P1-CR","P2-CL","P3-PR","P4-PL")
  write.table(pvals,paste("Datscan Results/",name,"Pvals",sep=""), sep=" ",
              quote = F, row.names = F, col.names = F)
}

#Calculate the correlation matrix
pea = cor(pvals[3:6],method = "pearson")  
if(MRI == 1){
write.table(pea,paste("MRI Results/WINSORIZE/",name,"/Cor",sep=""), sep=" ",
            quote = F, row.names = F, col.names = F)
} else {
  write.table(pea,paste("Datscan Results/",name,"Cor",sep=""), sep=" ",
              quote = F, row.names = F, col.names = F)
}
#----------------------------- PLOT THE RESULT OF TATES ---------------------------------------

if(MRI == 1){
res = read.table(paste("MRI Results/WINSORIZE/",name,"/Results_TATES.txt", sep=""))
} else {
res = read.table(paste("Datscan Results/",name,"Results_TATES.txt", sep=""))
}

# In order to plot the manhattan we need to add the bp column
whitebim = read.table("../Data/QC/White_2.bim",sep="")
BP = whitebim$V4
gwas = data.frame(res$V1,res$V2,BP,res$V4)
colnames(gwas) = c("CHR","SNP","BP","Pt")

library(qqman)
if(MRI == 1){
path = paste("MRI Results/WINSORIZE/",name,"/PLOT.png", sep="")
} else {
path = paste("Datscan Results/PLOT.png", sep="")
}

png(path, width = 24, height = 16, units = "in", res = 275)
par(mfrow=c(2,1))
manhattan(gwas, main = "Manhattan Plot of Pt - TATES", suggestiveline = -log10(1e-05), 
          genomewideline = -log10(5e-08), chr="CHR", bp="BP", p="Pt",ylim = c(0, 8), 
          col = c("red2", "green3") , annotatePval = 0.01)
qq(gwas$Pt, main = "Q-Q plot - TATES", col = "blue4")
dev.off()



