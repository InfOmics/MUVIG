info = read.csv("PPMIMERGE_20200419.csv")[-1]
id = read.csv("White_2_no_swedd.fam", sep=' ', header = F)

## from column 631 we have the imaging data infos

## taking only patients at baseline
library(tidyverse)
baseline = info %>% distinct(PATNO, .keep_all = TRUE)
df = baseline[baseline$PATNO %in% id[,1],]

lm1 <- lm(as.formula(paste(colnames(df)[632], "~", 
                           paste("APPRDX_enrol", "eTIV", "Age", "GENDER", sep ="+"),sep="")),
                            data = df)

min = c(summary(lm1)$coeff["APPRDX_enrol",4])
index = c(632)
names = colnames(df)[632]

for (i in 633:(ncol(df)-3)){
  
  z = df[,i]
  
  if (sum(z[!is.na(z)])>0){
  
    lm2 = lm(as.formula(paste(colnames(df)[i], "~", 
                            paste("APPRDX_enrol", "eTIV", "Age", "GENDER", sep ="+"),sep="")),
                            data = df)
    min = c(min,summary(lm2)$coeff["APPRDX_enrol",4])
    index = c(index,i)
    names = c(names, colnames(df)[i])
  }
}
features = data.frame(P = min, P.adj = p.adjust(min, method = "BH"), 
                      Index = index, Name = names)
features = features[order(features$P.adj, decreasing = F),]

write.table(features, "Ranked_Features.csv", quote = F, row.names = F, sep='\t')
