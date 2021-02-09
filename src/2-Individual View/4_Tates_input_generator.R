# This step is splitted in two parts, one using the MRI data and the other using the DATScan data

MRI = 1
if(MRI==1){
para_area_l = read.csv("../GWAS2/MRI/GWAS_MRI_lh_parahippocampal_area_Without_SWEDD_win.csv")
para_area_r = read.csv("../GWAS2/MRI/GWAS_MRI_rh_parahippocampal_area_Without_SWEDD_win.csv")
para_volume_r = read.csv("../GWAS2/MRI/GWAS_MRI_rh_parahippocampal_volume_Without_SWEDD_win.csv")
para_volume_l = read.csv("../GWAS2/MRI/GWAS_MRI_lh_parahippocampal_volume_Without_SWEDD_win.csv")
name = "Without SWEDD"
} else {
    caudate_r   = read.csv("../GWAS2/GWAS_DATScan_CAUDATE_R.csv")
    caudate_l   = read.csv("../GWAS2/GWAS_DATScan_CAUDATE_L.csv")
    putamen_r   = read.csv("../GWAS2/GWAS_DATScan_PUTAMEN_R.csv")
    putamen_l   = read.csv("../GWAS2/GWAS_DATScan_PUTAMEN_L.csv")
}

# Create the File of p-values
if(MRI==1){
# each P-value column correspond to a phenotype
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

# Calculate the correlation matrix
pea = cor(pvals[3:6],method = "pearson")  
if(MRI == 1){
write.table(pea,paste("MRI Results/WINSORIZE/",name,"/Cor",sep=""), sep=" ",
            quote = F, row.names = F, col.names = F)
} else {
  write.table(pea,paste("Datscan Results/",name,"Cor",sep=""), sep=" ",
              quote = F, row.names = F, col.names = F)
}
#----------------------------- PLOT THE RESULT OF TATES ------------------------------------

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

# to generate the plot you can just call the qqmna.r
source("qqman.r")
