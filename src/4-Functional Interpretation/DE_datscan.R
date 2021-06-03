#-------------------------------------------------------------------------------
#
# DE_putamen.R
# 
# Differential expression analysis on putamen dopamine uptake values.
#
# The below differential expression analysis test for genes showing a 
# differential expression by DaTscan putamen values, showed by the individuals
# participating to PPMI baseline RNA-seq.
# 
# The results are used as comparison metric for the other DE analysis, focusing 
# on searching differentially expressed genes by brain image values.
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

# plot DaTscan values distribution
p <- ggplot(metadata, aes(x = PUTAMEN_L)) +
  geom_density() + ggtitle("Putamen Left")
p
p <- ggplot(metadata, aes(x = PUTAMEN_R)) +
  geom_density() + ggtitle("Putamen Right")
p

# putamen normalization
metadata$PUTAMEN_L.norm <- RankNorm(metadata$PUTAMEN_L)
metadata$PUTAMEN_R.norm <- RankNorm(metadata$PUTAMEN_R)
p <- ggplot(metadata, aes(x = PUTAMEN_L.norm)) +
  geom_density() + ggtitle("Putamen Left Normalized")
p
p <- ggplot(metadata, aes(x = PUTAMEN_R.norm)) +
  geom_density() + ggtitle("Putamen Right Normalized")
p

####----------------------------------------------------------------------------
# PUTAMEN L

# Datscan values categorization
metadata$PUTAMEN_L.cat <- cut(
  metadata$PUTAMEN_L.norm, 
  breaks = c(quantile(metadata$PUTAMEN_L.norm, probs = seq(0, 1, by = 1/2))),
  labels = c(1,2)
)
metadata$PUTAMEN_L.cat[34] <- 1  # fix quantile() bug

##------------------------------------------------------------------------------
# DE analysis

# create DDS object
dds <- DESeqDataSetFromMatrix(
  countData = counts, colData = metadata, 
  design = ~ PUTAMEN_L.cat + gen + age_cat + ENROLL_CAT + educ
)

# filter by expression
keep <- apply(counts(dds), 1, function(x){sum(x != 0) > 132})
dds <- dds[keep,]

# normalize counts
vsd <- vst(dds, blind = F)
p <- plotPCA(vsd, intgroup = "gen", returnData=F)
p + geom_label(aes(label = colnames(dds))) + 
  ggtitle("Gender")  #  batch effect by gender

# correct counts to remove batch effect
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$gen)
p <- plotPCA(vsd, intgroup = "gen", returnData=F)
p + geom_label(aes(label = colnames(dds))) + 
  ggtitle("Gender")  # batch effect removed

# plot variability via PCA
p <- plotPCA(vsd, intgroup = "ENROLL_CAT", returnData=F)
p + geom_label(aes(label = colnames(dds))) + 
  ggtitle("Enrollment category")

p <- plotPCA(vsd, intgroup = "PUTAMEN_L.cat", returnData=F)
p + geom_label(aes(label = colnames(dds))) + 
  ggtitle("Putamen Left value category")

# DESeq2 analysis
dds <- DESeq(dds)

# compute results
results.putamen.l <- results(dds, contrast = c("PUTAMEN_L.cat", "1", "2"))
results.putamen.l <- results.putamen.l[!is.na(results.putamen.l$padj),]
results.putamen.l <- results.putamen.l[results.putamen.l$padj < .1,]
dim(results.putamen.l)


##------------------------------------------------------------------------------
# enrichment analysis

# map ENSEMBL IDs to ENTREZ IDs
ensembl.ids <- rownames(results.putamen.l)
ensembl.ids <- sapply(
  rownames(results.putamen.l),
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
results.putamen.l <- results[entrez.ids$ENSEMBL,]  # remove unmapped genes
rownames(results.putamen.l) <- entrez.ids$ENTREZID

# enrichment with fgsea - reactome db
lfc <- results.putamen.l$log2FoldChange
names(lfc) <- rownames(results.putamen.l)
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
  "../data/RNA-seq-data/pathways-putamen-l-reactome.tsv",
  sep = "\t",
  sep2 = c("", " ", "")
)

# enrichment with fgsea - kegg db
pathways <- gmtPathways("../data/RNA-seq-data/GMT/c2.cp.kegg.v7.4.entrez.gmt")
enr.fgsea.kegg <- fgsea(
  pathways, 
  lfc,
  nPermSimple = 10000,
  eps = 0,
  minSize = 10
)

# subset to signifcant pathways (q-value < 0.1)
enr.fgsea.kegg <- enr.fgsea.kegg[enr.fgsea.kegg$padj < .1,]

# store results
fwrite(
  as.data.frame(enr.fgsea.kegg),
  "../data/RNA-seq-data/pathways-putamen-l-kegg.tsv",
  sep = "\t",
  sep2 = c("", " ", "")
)



