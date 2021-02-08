
datscan = read.csv("../DATScan_Analysis.csv")
#Common in statstical genetics is to use rank-based inverse normal transformation
library(RNOmni)
putamen_norm_r = rankNorm(datscan$PUTAMEN_R) #rank based inverse normal transformation (INT) 
putamen_norm_l = rankNorm(datscan$PUTAMEN_L)
caudate_norm_r = rankNorm(datscan$CAUDATE_R)
caudate_norm_l = rankNorm(datscan$CAUDATE_L)

# basta fare apply(datscan,2,rankNorm)

datascan_norm = data.frame(datscan$PATNO,datscan$EVENT_ID,datscan$SCAN_DATE,caudate_norm_r
                           ,caudate_norm_l,putamen_norm_r,putamen_norm_l)
colnames(datascan_norm) <- c("PATNO","EVENT_ID","SCAN_DATE","CAUDATE_R","CAUDATE_L",
                             "PUTAMEN_R","PUTAMEN_L")
write.csv(datascan_norm,"../DATScan_Norm_complete.csv",quote=F,row.names = F)
