##############################################################################################
#
# 2-RNASeq.R
#
# The following R script compute a Differential expression analysis using the count data
# as input, performed using LIMMA and variancePartition. LIMMA builds a linear model to fit 
# gene expression-phenotype association. At the end will be extract a table of the top-ranked 
# genes from a linear model fit. To increase the time execution of the algorithm the core
# part of the method is executed using multiple core with the SnowParam function of the 
# BiocParallel package. 
#
##############################################################################################


# R packages required to run the script:
# - edgeR
# - variancePartition
# - BiocParallel
# - tidyverse
# - AnnotationDbi
# - AnnotationHub
# - RColorBrewer
  

if(!require("edgeR", character.only = TRUE))
{
  BiocManager::install("edgeR")
  if(!require("edgeR", character.only = TRUE))
  {
    stop("edgeR package not found")
  }
}

if(!require("variancePartition", character.only = TRUE))
{
  BiocManager::install("variancePartition")
  if(!require("variancePartition", character.only = TRUE))
  {
    stop("variancePartition package not found")
  }
}

if(!require("BiocParallel", character.only = TRUE))
{
  BiocManager::install("BiocParallel")
  if(!require("BiocParallel", character.only = TRUE))
  {
    stop("BiocParallel package not found")
  }
}


if(!require("tidyverse", character.only = TRUE))
{
  install.packages("tidyverse")
  if(!require("tidyverse", character.only = TRUE))
  {
    stop("tidyverse package not found")
  }
}

if(!require("AnnotationDbi", character.only = TRUE))
{
  BiocManager::install("AnnotationDbi")
  if(!require("AnnotationDbi", character.only = TRUE))
  {
    stop("AnnotationDbi package not found")
  }
}

if(!require("AnnotationHub", character.only = TRUE))
{
  BiocManager::install("AnnotationHub")
  if(!require("AnnotationHub", character.only = TRUE))
  {
    stop("AnnotationHub package not found")
  }
}  
 
if(!require("RColorBrewer", character.only = TRUE))
{
  install.packages("RColorBrewer")
  if(!require("RColorBrewer", character.only = TRUE))
  {
    stop("RColorBrewer package not found")
  }
}    
  
# load the required packages
suppressPackageStartupMessages(c( 
                                  library(edgeR),
                                  library(variancePartition),
                                  library(BiocParallel),
                                  library(tidyverse),
                                  library(AnnotationDbi),
                                  library(AnnotationHub),
                                  library(RColorBrewer),
                                  require('gtools') ))



# First of all we need to create a clean dataframe i.e. without NA 
# and containing all the imaging and confounding information to be used

# list of counts data
files = list.files("Final_422/")
files = mixedsort(files) # order the count files

# load the different information
cov     = read.csv("covariate_mri.txt", sep=" ")[-c(2:12)]
datscan = read.csv("phenotype_datscan.txt", sep=" ")[-2]
mri     = read.csv("phenotype_mri.txt", sep=" ")[-c(2,3)]
sex     = read.csv("PPMI_eu_noswedd_ds.fam", sep=" ", header=F)[-c(2:4,6)] 
colnames(sex) = c("FID","SEX")

# generate a single dataframe containing all the information needed
# to run the following line is necessary to have the same column name for all the dataframe
info = Reduce(merge, list(cov,datscan,mri,sex))
# delete the NA values 
info = na.omit(info)

# Since we now have 388 subjects we need to select the count data only for those subjects
PatientsIDs = files[ which(
                      sapply( strsplit(files, ".", fixed = TRUE), "[[", 2) 
                        %in% info$FID == TRUE )  ]


# read the data into a DGEList object
x <- readDGE(PatientsIDs, "Final_422/", columns=c(1,7), header = T, skip=1)
# can be useful to store the DGEList object to a local file in order to reload it in a faster way
saveRDS(x, "DGEList.rds")

x = readRDS("DGEList.rds")

# see the dimensions
dim(x)

# Since x now have multiple row per patient we need to extend in the same way the dataframe
# containing the information

multipleIds = as.data.frame(
                  sapply( 
                      strsplit(files, ".", fixed = TRUE), "[[", 2) )

colnames(multipleIds) <- "FID"
info = merge(info, multipleIds, by="FID")

# adding all the information into the sample of the DGEList object
x$samples[,5:18] = info[,c(6:15,1:3,16)]
x$samples$group  = as.factor( ifelse(info$ENROLL_CAT == 1 , "HC", "PD") )

######################## 1. ORGANISING GENE ANNOTATIONS ############################ 
  
geneid <- rownames(x) 
geneid <- gsub("\\..*","", geneid) # removing the dot (version) after the gene name

# convert the gene ensemble into entrezid and symbol
ah <- AnnotationHub()
# retrieve the latest version of EnsDb
ens   = query(ah, c("EnsDb", "Hsapiens", "103"))[[1]]
genes = select(ens, key=geneid, columns=c("ENTREZID", "SYMBOL"), keytype="GENEID")
# As with any gene ID, Entrez gene IDs may not map one-to-one to the gene information of interest. 
# It is important to check for duplicated gene IDs
genes <- genes[!duplicated(genes$GENEID),] # removing duplicates
x$genes <- genes # adding the gene name to the DGEList object since now we have fewer genes 


######################## 2. DATA PRE-PROCESSING ############################

# 2.1 TRANSFORMATION FROM THE RAW-SCALE

  # It is common practice to transform raw counts onto a scale that accounts for 
  # library size differences 
  
  # since we do not have the gene length information we use the counts per million (CPM) 
  # transformation
  cpm = cpm(x)

# 2.2 REMOVING GENES THAT ARE LOWLY EXPRESSED

  # Some genes are unexpressed throughout all samples
  t = table(rowSums(x$counts==0)==ncol(x))
  zeroCountsPer = (t[2] / (t[1] + t[2]) ) * 100
  # We are lucky because only the 1.6% of genes in this dataset have zero counts across all samples

  # Determine which genes have sufficiently large counts to be retained in a statistical analysis
  keep.exprs <- rowSums(cpm(x)>0.4) >= 5
  
  keep.exprs <- filterByExpr(x, group = x$samples$group)
  x_filtered <- x[keep.exprs,]
  
  # create some useful plot before and after the filtering
  # Plotting the distribution log-CPM values shows that a sizeable proportion of genes 
  # within each sample are either unexpressed or lowly-expressed with log-CPM values that 
  # are small or negative 
  L <- mean(x$samples$lib.size)   * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  lcpm.cutoff <- log2(10/M + 2/L)
  
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  
  par(mfrow=c(1,2))
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
  lcpm <- cpm(x_filtered, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")
  abline(v=lcpm.cutoff, lty=3)
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }

# 2.3 NORMALISING GENE EXPRESSION DISTRIBUTIONS

  # Normalize: Normalisation is required to ensure that the expression distributions
  # of each sample are similar across the entire experiment.
  
  # boxplot, is useful in determining whether any samples are dissimilar to others
  #To give a better visual representation of the effects of normalisation, the data was 
  # duplicated then adjusted so that the counts of the first sample are reduced to 5% of 
  # their original values, and in the second sample they are inflated to be 5-times larger.
  x2 = x_filtered
  x2$samples$norm.factors = 1
  x2$counts[,1] = ceiling(x2$counts[,1]*0.05)
  x2$counts[,2] = x2$counts[,2]*5
  
  lcpm_unnorm = cpm(x2, log=TRUE)
  tmp   = calcNormFactors(x2, method = "TMM")  #instead of TMM you can use different method
  lcpm_norm  = cpm(tmp, log=TRUE)
  
  boxplot(lcpm_unnorm[,1:50], las=2, col=col, main="")
  title(main="A. Example: Unnormalised data",ylab="Log-cpm")
  
  boxplot(lcpm_norm[,1:40], las=2, col=col, main="")
  title(main="B. Example: Normalised data - TMM",ylab="Log-cpm")

  # normalization, through the various plots the TMMwsp method appears to be the best
  x_norm <- calcNormFactors(x_filtered, method = "TMMwsp")

######################### 3 DIFFERENTIAL EXPRESSION ANALYSIS ########################### 

param = SnowParam(4, "FORK", progressbar=TRUE) # run in parallel using 4 core
register(param)

for (i in colnames(x_norm$samples)[c(5:14)]){
  factors = c("~0",i,"eTIV","AGE","SEX","(1|FID)")
  form = as.formula(paste(factors, collapse="+"))
  
  # Removing heteroscedascity from count data
  # voomWithDreamWeights() replaces voom() to estimate precision weights
  vobjDream = voomWithDreamWeights(x_norm, form, x_norm$samples, save.plot = T)
  
  # save the Mean Variance Trend plot
  png(paste0("PLOT/MeanVarianceTrend_",i,".png"),width = 16, height = 12, units="in", res=100) 
  plot(vobjDream$voom.xy, main ="Voom: Mean Variance Trend", 
       ylab="Sqrt (standard deviation )", xlab = "log2( count size + 0.5 )")
  lines(vobjDream$voom.line, col="red")
  dev.off() 
  
  # dream() replaces lmFit() to estimate regression coefficients.
  fitmm = dream(vobjDream, form, x_norm$samples, BPPARAM=param)
  
  # The contrast of the linear model of interest is the phenotype
  tt <- topTable(fitmm, n = Inf, coef=i , sort.by = "P", adjust.method="BH") 
  path = paste0("Results/Genes_", colnames(fitmm)[1], ".csv")
  write.csv(tt, file=path, quote = F, row.names = F)
  
  print(paste("Table created for the feature: ", i))
}




