library(RNOmni)
library(edgeR)
library(variancePartition)
library(VennDiagram)
library(BiocParallel)
library(fgsea)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(Homo.sapiens)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v75)
require('gtools')

# create the dataframe with all the infos
files_no_swedd = list.files("./Data/Count/With_SWEDD/")
files_no_swedd = mixedsort(files_no_swedd)
cov = read.csv("Data/Info/With_Swedd/Covariate_With_SWEDD.txt", sep="\t")
datscan = read.csv("Data/Info/With_Swedd/Phenotype_DAT.txt", sep="\t")
mri = read.csv("Data/Info/With_Swedd/Phenotype_With_SWEDD_win_MRI.txt", sep="\t")
fam = read.csv("Data/Info/With_Swedd/White_2.fam", sep=" ", header=F)[-c(2:4,6)]
colnames(fam) = c("FID","SEX")

# retrieve a list of all id repeated for each patients
my.list <- strsplit(files_no_swedd, ".", fixed = TRUE)
IDs <- sapply(my.list, "[[", 2)

tmp = merge(cov,datscan,by="FID")
info = merge(tmp,mri,by="FID")
columnsOut = which(colnames(info) %in% c("IID","IID.y","APPRDX2.y"))
info = info[,-columnsOut]
colnames(info)[c(2,10)] = c("IID","GROUP")
info = merge(info,fam,by="FID")

tmp2 = as.data.frame(IDs)
colnames(tmp2) <- "FID"
info_multiple = merge(info, tmp2, by="FID")

# delete the NA values 
delete = rownames(info_multiple)[!complete.cases(info_multiple$eTIV)]
info_filtered = info_multiple[!rownames(info_multiple) %in% delete,]

lista = IDs %in% info_filtered$FID
indici = which(lista==TRUE)
nostriPazienti = files_no_swedd[indici]

x <- readDGE(nostriPazienti, "./Data/Count/With_SWEDD/", columns=c(1,7), header = T, skip=1)
dim(x)

x$samples$lh_parahippocampal_area   <- info_filtered$lh_parahippocampal_area
x$samples$rh_parahippocampal_area   <- info_filtered$rh_parahippocampal_area
x$samples$lh_parahippocampal_volume <- info_filtered$lh_parahippocampal_volume
x$samples$rh_parahippocampal_volume <- info_filtered$rh_parahippocampal_volume
x$samples$putamen_r  <- info_filtered$PUTAMEN_R
x$samples$putamen_l  <- info_filtered$PUTAMEN_L
x$samples$caudate_r  <- info_filtered$CAUDATE_R
x$samples$caudate_l  <- info_filtered$CAUDATE_L
x$samples$gender     <- info_filtered$SEX
x$samples$eTIV       <- info_filtered$eTIV
x$samples$age        <- info_filtered$AGE
x$samples$individual <- info_filtered$FID # use later to better weight the multiple counts

# Organising gene annotations
geneid <- rownames(x) 
geneid <- gsub("\\..*","", geneid) # removing the dot (version) after the gene name

genes = select(EnsDb.Hsapiens.v75, key=geneid, columns=c("ENTREZID", "SYMBOL"), keytype="GENEID")
genes <- genes[!duplicated(genes$GENEID),] # removing duplicates
x$genes <- genes

# Determine which genes have sufficiently large counts to be retained in a statistical analysis
keep.exprs <- rowSums(cpm(x)>1) >= 5
x <- x[keep.exprs,] 

# Normalize: Normalisation is required to ensure that the expression distributions
# of each sample are similar across the entire experiment.
x <- calcNormFactors(x, method = "TMM")


# ---------------------Differential expression analysis-------------------------------

param = SnowParam(6, "FORK", progressbar=TRUE) # run in parallel using 6 core
register(param)

for (i in colnames(x$samples)[c(5:12)]){
factors = c("~0",i,"eTIV","age","gender","(1|individual)")
form = as.formula(paste(factors, collapse="+"))
#form <- ~0 + caudate_r + eTIV + age + gender + (1|individual)
vobjDream = voomWithDreamWeights(x, form, x$samples, plot=F)
print(paste("vobjDream object created for:",i))
fitmm = dream(vobjDream, form, x$samples, BPPARAM=param)
print(paste("fitmm object created for:",i))
tt <- topTable(fitmm, n = Inf, coef=i , sort.by = "P", adjust.method="BH") 
# contrast of the linear model is of interest = 1
path = paste0("Results/With_Swedd/TopTable_", colnames(fitmm)[1], ".csv",sep="")
write.csv(tt, file=path, quote = F, row.names = F)
print(paste("Table created for:",i))
}

################################## Zip a file ################################## 
files2zip <- dir('/home/ubuntu/Results/With_Swedd/', full.names = TRUE)
zip(zipfile = 'Results_With_Swedd', files = files2zip)
################################################################################


