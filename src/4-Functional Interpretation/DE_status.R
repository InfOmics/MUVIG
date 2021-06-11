#-------------------------------------------------------------------------------
#
# DE_status.R
# 
# Differential expression analysis on subjects status
#
# The below differential expression analysis test for genes showing a 
# differential expression by subject status.
# 
# The results are used as comparison metric for the other DE analysis, focusing 
# on searching differentially expressed genes by brain image values.
#
# After the differential expression analysis we'll be carried out en enrichment
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
#   - data.table

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

if(!require("forg.Hs.eg.db", character.only = TRUE))
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

# retrieve patno from counts file names
getid.master <- function(s){
  patno <- unlist(strsplit(s, ".", fixed = T))[4]
  return(patno)
}
samples.ids <- colnames(counts)
samples.ids.master <- sapply(samples.ids, getid.master)
names(samples.ids.master) <- NULL
colnames(counts) <- samples.ids.master

# sort metadata
rownames(metadata) <- metadata$PATNO
metadata <- metadata[colnames(counts),]

# cast vars to factors for DE analysis
metadata$gen <- as.factor(metadata$gen)
metadata$age_cat <- as.factor(metadata$age_cat)
metadata$educ <- as.factor(metadata$educ)


##------------------------------------------------------------------------------
# DE analysis

# create DDS object from counts matrix
dds <- DESeqDataSetFromMatrix(
  countData = counts, 
  colData = metadata, 
  design = ~ gen + age_cat + educ + ENROLL_CAT
)

# filter by expression
keep <- apply(counts(dds), 1, function(x){sum(x != 0) > 175})
dds <- dds[keep,]

# normalize counts
vsd <- vst(dds, blind = F)
p <- plotPCA(vsd, intgroup = "gen", returnData = F)
p + geom_label(aes(label = colnames(dds))) +
  ggtitle("Gender")  # sex creates considerable batch effect

# correct to remove sex batch effect
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$gen) 
p <- plotPCA(vsd, intgroup = "gen", returnData = F)
p + geom_label(aes(label = colnames(dds))) +
  ggtitle("Gender")

# plot variability via PCA
p <- plotPCA(vsd, intgroup = "ENROLL_CAT", returnData = F)
p + geom_label(aes(label = colnames(dds))) +
  ggtitle("Enrollment category")

# de analysis
dds <- DESeq(dds, blind = F)

# compute results
results <- results(dds, contrast = c("ENROLL_CAT", "HC", "PD"))
results <- results[!is.na(results$padj),]
results <- results[results$padj < .1,]
dim(results)  # 2201 de genes


##------------------------------------------------------------------------------
# enrichment analysis

# map ENSEMBL IDs to ENTREZ IDs
ensembl.ids <- rownames(results)
ensembl.ids <- sapply(
  rownames(results),
  function(x){
    unlist(strsplit(x, fixed = T, split = "."))[1]
  }
)
names(ensembl.ids) <- NULL  # remove old names 
entrez.ids <- select(
  org.Hs.eg.db, 
  keys = ensembl.ids, 
  columns = "ENTREZID",
  keytype = "ENSEMBL"
)
entrez.ids <- entrez.ids[!duplicated(entrez.ids$ENSEMBL),]  # remove duplicates
entrez.ids <- entrez.ids[!is.na(entrez.ids$ENTREZID),]  # remove NAs
results <- results[entrez.ids$ENSEMBL,]  # remove unmapped genes
rownames(results) <- entrez.ids$ENTREZID

# enrichment with fgsea - reactome db
lfc <- results$log2FoldChange
names(lfc) <- rownames(results)
pathways <- reactomePathways(names(lfc))  # recover pathways with de genes
enr.fgsea.reactome <- fgsea(
  pathways, 
  lfc,
  nPermSimple = 10000,
  eps = 0,
  minSize = 10
)

# subset to signifcant pathways (q-value < 0.1)
enr.fgsea.reactome <- enr.fgsea.reactome[enr.fgsea.reactome$padj < .1,]

# store results
fwrite(
  as.data.frame(enr.fgsea.reactome),
  "../data/RNA-seq-data/pathways-cond-reactome.tsv",
  sep = "\t",
  sep2 = c("", " ", "")
)

## -- add nbe -- ##
