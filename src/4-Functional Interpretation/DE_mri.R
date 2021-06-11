#-------------------------------------------------------------------------------
#
# DE_mri.R
# 
# Differential expression analysis on parahippocampal morphoogical measures 
# recovered from brain MRI images. The analyzed features are parahippocampl
# volume, area and thickness.
#
# The below differential expression analysis test for genes showing a 
# differential expression by MRI parahippocampal morphological values. 
# 
#
# After the differential expression analysis we'll be carried out an enrichment
# analysis to recover the  biological pathways potentially perturbed by 
# differentially expressed genes.
#


# R packages required to run the script:
#   - DESeq2
#   - ggplot2
#   - limma
#   - fgsea
#   - org.Hs.eg.db
#   - ReactomePA


if(!require("BiocManager", character.only = TRUE))
{
  install.packages("BiocManager")
  if(!require("BiocManager", character.only = TRUE))
  {
    stop("BiocManager package not found")
  }
}

if(!require("DESeq2", character.only = TRUE))
{
  BiocManager::install("DESeq2")
  if(!require("DESeq2", character.only = TRUE))
  {
    stop("DESeq2 package not found")
  }
}

if(!require("ggplot2", character.only = TRUE))
{
  install.packages("ggplot2")
  if(!require("ggplot2", character.only = TRUE))
  {
    stop("ggplot2 package not found")
  }
}

if(!require("limma", character.only = TRUE))
{
  BiocManager::install("limma")
  if(!require("limma", character.only = TRUE))
  {
    stop("limma package not found")
  }
}

if(!require("fgsea", character.only = TRUE))
{
  BiocManager::install("fgsea")
  if(!require("fgsea", character.only = TRUE))
  {
    stop("fgsea package not found")
  }
}

if(!require("org.Hs.eg.db", character.only = TRUE))
{
  BiocManager::install("org.Hs.eg.db")
  if(!require("org.Hs.eg.db", character.only = TRUE))
  {
    stop("org.Hs.eg.db")
  }
}

if(!require("ReactomePA", character.only = TRUE))
{
  BiocManager::install("ReactomePA")
  if(!require("ReactomePA", character.only = TRUE))
  {
    stop("ReactomePA package not found")
  }
}

if(!require("data.table", character.only = TRUE))
{
  install.packages("data.table")
  if(!require("data.table", character.only = TRUE))
  {
    stop("data.table package not found")
  }
}


# load required packages
library(DESeq2)
library(ggplot2)
library(limma)
library(fgsea)
library(org.Hs.eg.db)
library(ReactomePA)
library(data.table)


##------------------------------------------------------------------------------
# load data

counts <- read.table("../data/RNA-seq-data/PPMI-RNA-seq_baseline_rawCounts.txt")
metadata <- read.csv("../data/RNA-seq-data/baseline_metadata.csv")

# retrieve patno from counts matrix
getid <- function(x){
  patno <- unlist(strsplit(s, ".", fixed = T))[4]
  return(patno)
}
samples.ids <- colnames(counts)
sample.ids.new <- sapply(samples.ids, getid)
names(sample.ids.new) <- NULL  # remove old subjects names
colnames(counts) <- sample.ids.new
View(counts)

# filter metadata
rownames(metadata) <- metadata$PATNO
metadata <- metadata[colnames(counts),]  # keep data related to considered subjects
View(metadata)

# cast useful vars to factor
metadata$gen <- as.factor(metadata$gen)  # gender
metadata$age_cat <- as.factor(metadata$age_cat)  # age category
metadata$educ <- as.factor(metadata$educ)  # education years

# explore and normalize intracranial area (eTIV)
p <- ggplot(metadata, aes(x = eTIV)) +
  geom_density() + ggtitle("eTIV")
p
metadata$eTIV.norm <- RankNorm(metadata$eTIV)
p <- ggplot(metadata, aes(x = eTIV.norm)) + 
  geom_density() + ggtitle("eTIV normalized")

# categorization of eTIV in three levels
metadata$eTIV.cat <- cut(
  metadata$eTIV.norm,
  breaks = c(quantile(metadata$eTIV.norm, probs = seq(0, 1, by = 1/3))),
  labels = c(1,2,3)
)
metadata[71] <- 1

# parahippocampal area normalization
p <- ggplot(metadata, aes(x = lh_parahippocampal_area)) +
  geom_density() + ggtitle("Parahippocampal area left")
p
metadata$lh_parahippocampal_area.norm <- RankNorm(
  metadata$lh_parahippocampal_area
)
p <- ggplot(metadata, aes(x = lh_parahippocampal_area.norm)) +
  geom_density() + ggtitle("Normalized parahippocampal area left")
p

####----------------------------------------------------------------------------
# PARAHIPPOCAMPAL AREA L

print("Performing DE analysis on PARAHIPPOCAMPAL Area Left", quote = F)

# Parahippocampal area left values categorization
metadata$lh_parahippocampal_area.cat <- cut(
  metadata$lh_parahippocampal_area.norm, 
  breaks = c(
    quantile(metadata$lh_parahippocampal_area.norm, probs = seq(0, 1, by = 1/2))
  ),
  labels = c(1,2)
)
metadata$lh_parahippocampal_area.cat[107] <- 1  # fix quantile() bug

##------------------------------------------------------------------------------
# DE analysis

dds <- DESeqDataSetFromMatrix(
  countData = counts, colData = metadata,
  design = ~ lh_parahippocampal_area.cat + eTIV.cat + gen + age_cat + ENROLL_CAT + educ
)

# filter by expression
keep <- apply(counts(dds), 1, function(x){sum(x != 0) > 133})
dds <- dds[keep,]

# normalize counts
vsd <- vst(dds, blind = F)
p <- plotPCA(vsd, intgroup = "gen", returnData = F)
p + geom_label(aes(label = colnames(dds))) +
  gggtitle("Gender")  # BATCH effect
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$gen)
p <- plotPCA(vsd, intgroup = "gen", returnData = F)
p + geom_label(aes(label = colnames(dds))) +
  ggtitle("Gender")  # batch effect removed

# plot variability via PCA
p <- plotPCA(vsd, intgroup = "ENROLL_CAT", returnData = F)
p + geom_label(aes(label = colnames(dds))) + 
  gggtitle("Enrollment category")
p <- plotPCA(vsd, intgroup = "eTIV.cat", erturnData = F)
p + geom_label(aes(label = colnames(dds))) + 
  ggtitle("eTIV category")
p <- plotPCA(vsd, intgroup = "lh_parahippocampal_area.cat", returnData = F)
p + geom_label(aes(label = colnames(dds))) + 
  ggtitle("Parahippocampal area left category")

# DESeq2 analysis

dds <- DESeq(dds)

# contrast by parahippocampal left area category
results.pharea.l <- results(
  dds, 
  contrast = c("lh_parahippocampal_area.cat", "1", "2")
)
results.pharea.l <- results.pharea.l[!is.na(results.pharea.l$p.adj),]
results.pharea.l <- results.pharea.l[results.pharea.l$p.adj < .1,]
dim(results.pharea.l)




