library(devtools)
options(stringsAsFactors = FALSE) 
library(grace) 
library(fda) 
library(ggplot2) 
library(mvtnorm) 
library(tidyverse)
library(qqman)

data = read.csv("../../Processed_Images/PPMIMERGE_20200419.csv")[-1]
# We take only the features of interest
data = data[,c(1,2,5,645,679,783,818)]
data = data[complete.cases(data),]
# We select just the first row for each patient
data = data %>% distinct(PATNO, .keep_all = TRUE)

datscan = read.csv("../../Imaging/DATScan_Norm.csv")
datscan = datscan[,c(1,4:ncol(datscan))]
datscan = datscan %>% distinct(PATNO, .keep_all = TRUE)

pheno = read.csv("../GWAS2/MRI/Phenotype_Without_SWEDD_win.txt",sep="\t")
cov =   read.csv("../GWAS2/MRI/Covariate_Without_SWEDD.txt",sep="\t")
pheno = pheno[complete.cases(pheno),]
cov = cov[complete.cases(cov),]

# sort of normalization for the value in the phenotype
for (i in 4:ncol(pheno)){
  x = sd(pheno[[i]])
  pheno[i] = (pheno[[i]] - min(pheno[[i]])) / x
}

#  ****************************** USING MRI DATA ******************************
df = data.frame("PATNO"=pheno$FID,
                "AGE"=cov$AGE, # to simulate a follow up data 
                "GROUP"=pheno$APPRDX2,
                "Para_R_Area"=pheno$rh_parahippocampal_area,
                "Para_L_Area"=pheno$lh_parahippocampal_area,
                "Para_R_Volume"=pheno$rh_parahippocampal_volume,
                "Para_L_Volume"=pheno$lh_parahippocampal_volume)
z = c("Para_R_Area","Para_L_Area","Para_R_Volume","Para_L_Volume")
tmp = df[,4:ncol(df)]
df_grace = data.frame("id" = rep(df$PATNO,each=4),
                      "argvals" = rep(df$AGE,each=4),
                      "group" = rep(df$GROUP,each=4),
                      "Outcome" = rep(z,nrow(df)),
                      "Y" = as.numeric(as.list(t(tmp))))

grace.simulation.fits <- grace(id = df_grace$id,
                               argvals = df_grace$argvals,
                               y = df_grace$Y,
                               outcome = df_grace$Outcome,
                               group = df_grace$group,
                               plots = T)

ggplot(grace.simulation.fits$sigma, aes(Iteration, sigma, group=outcome))+
  geom_line(aes(colour = outcome)) +
  ylab(expression(sigma[epsilon]))

for (i in z){
  tmp <- grace.simulation.fits$fits[[i]]$subset
  write.csv(tmp, paste("MRI/WINSORIZE/NO SWEDD/Grace_Results_",i,".csv",sep=""),row.names = F,quote=F)
}
  
#  ****************************** USING DATSCAN DATA  ******************************

df2 = data.frame("PATNO"=swedd$PATNO,
                 "AGE"=swedd$Age,
                 "GROUP"=swedd$APPRDX_enrol,
                 "Caudate_R"=swedd$CAUDATE_R,
                 "Caudate_L"=swedd$CAUDATE_L,
                 "Putamen_R"=swedd$PUTAMEN_R,
                 "Putamen_L"=swedd$PUTAMEN_L)

v = c("Caudate_R","Caudate_L","Putamen_R","Putamen_L")

for (i in v){
  df2[i] = (df2[[i]] - min(df2[[i]],na.rm=T))/sd(df2[[i]],na.rm = T)
}

tmp2 = df2[,4:ncol(df2)]
df_grace2 = data.frame("id" = rep(df2$PATNO,each=4),
                       "argvals" = rep(df2$AGE,each=4),
                       "group" = as.factor(rep((ifelse(df2$GROUP==2,"HC","PD")),each=4)),
                       "Outcome" = rep(v,395),
                       "Y" = as.numeric(as.list(t(tmp2))))

grace.fits <- grace(id = df_grace2$id,
                    argvals = df_grace2$argvals,
                    y = df_grace2$Y,
                    outcome = df_grace2$Outcome,
                    group = df_grace2$group,
                    plots = T)

ggplot(grace.fits$sigma, aes(Iteration, sigma, group=outcome))+
  geom_line(aes(colour = outcome)) +
  ylab(expression(sigma[epsilon]))

for (i in v){
  tmp <- grace.fits$fits[[i]]$subset
  write.csv(tmp, paste("DATSCAN/WITHOUT SWEDD/Grace_Results_",i,".csv",sep=""),row.names = F,quote=F)
}

#  ****************************** GWAS ANALYSIS ******************************

# conduct a gwas analysis using as phenotype the gamma values returned by grace instead of 
# using the standard P-values

# path to the executable
path = "/Tools/plink_mac_20191028/plink"
                  
df = read.csv("Grace_Results_Para_L_Area.csv")
df = df[df$id %in% pheno$FID,]

pheno = data.frame(df$id, cov$IID, df$group, df$gamma0)
colnames(pheno) <- c("FID","IID","GROUP","GAMMA")
pheno[is.na(pheno)] <- "NA"
cov[is.na(cov)] <- "NA"
write.table(pheno, "pheno.csv", quote = F, row.names = F, sep='\t')
write.table(cov, "cov.csv", quote = F, row.names = F, sep='\t')

name = "MRI"
cmd <- sprintf(" --bfile ../Data/QC/White_2 --pheno pheno.csv --covar cov.csv --pheno-name GAMMA --covar-name AGE,PC1-PC5 --linear sex --adjust --out %s" 
               , name)
system(paste0(path,cmd))

gwas = read.table(sprintf("%s.assoc.linear", name), header = T, sep="")
gwas = gwas[gwas$TEST=="ADD",]
write.table(gwas, sprintf("GWAS_%s.csv", name), quote = F, row.names = F)

png(sprintf("%s_Manhattan_QQ_NO_SWEDD_win.png", name), width = 24, height = 16, units = "in", res = 275)
par(mfrow=c(2,1)) 
manhattan(gwas, main = sprintf("Manhattan Plot of %s", name), suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), 
          chr="CHR", bp="BP", p="P",ylim = c(0, 8), col = c("red2", "green3") , annotatePval = 0.01)
qq(gwas$P, main = sprintf("Q-Q plot of %s", name), col = "blue4")
dev.off()
