####################################################################################
#
# 1-Annotation.R
#
# The following R script generate a table containing different information
# related to the resulting significant SNPs obtained after applying single-trait
# based statistical significance of GWAS results with TATES and the growth curve
# models computed using GRACE. This method uses NCBI's refSNP for information 
# related to the latest dbSNP build and the STRING database for known and
# predicted protein-protein interaction. The interactions include direct (physical) 
# and indirect (functional) associations. It also adds the pathways's list involving 
# a gene, using data from WikiPathways.
# 
# Information include: 
#   - query: The rs ID that was queried
#   - chromosome
#   - bp
#   - class: The rsid's 'class' (SNP class description here: 
#            https://www.ncbi.nlm.nih.gov/projects/SNP/snp_legend.cgi?legend=snpClass)
#   - rsid: Reference SNP cluster ID
#   - gene (SYMBOL): NA if the SNP does not map in a gene
#   - alleles
#   - ancestral_allele: allele as described in the current assembly
#   - variation_allele: difference to the current assembly
#   - seqname: chromosome RefSeq reference
#   - assembly
#   - hgvs
#   - ref_seq
#   - minor: The allele for which the MAF is computed
#   - maf:   The minor allele frequency of the SNP
#   - PPI:   The Protein-Protein interaction network list for the gene, 
#            NA if the SNP does not map in a gene
#   - Pathway: The Pathway list for the gene, 
#              NA if the gene can not be mapped
#
####################################################################################


# R packages required to run the script:
# - rsnps
# - httr
# - BiocManager
# - clusterProfiler
# - org.Hs.eg.db
# - rWikiPathways

if(!require("rsnps", character.only = TRUE))
{
  install.packages("rsnps")
  if(!require("rsnps", character.only = TRUE))
  {
    stop("rsnps package not found")
  }
}

if(!require("httr", character.only = TRUE))
{
  install.packages("httr")
  if(!require("httr", character.only = TRUE))
  {
    stop("httr package not found")
  }
}

if(!require("BiocManager", character.only = TRUE))
{
  install.packages("BiocManager")
  if(!require("BiocManager", character.only = TRUE))
  {
    stop("BiocManager package not found")
  }
}

if(!require("clusterProfiler", character.only = TRUE))
{
  BiocManager::install("clusterProfiler")
  if(!require("clusterProfiler", character.only = TRUE))
  {
    stop("clusterProfiler package not found")
  }
}

if(!require("org.Hs.eg.db", character.only = TRUE))
{
  BiocManager::install("org.Hs.eg.db")
  if(!require("org.Hs.eg.db", character.only = TRUE))
  {
    stop("org.Hs.eg.db package not found")
  }
}

if(!require("rWikiPathways", character.only = TRUE))
{
  BiocManager::install("rWikiPathways")
  if(!require("rWikiPathways", character.only = TRUE))
  {
    stop("rWikiPathways package not found")
  }
}

# load the required packages
suppressPackageStartupMessages(library(rsnps))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(BiocManager))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(rWikiPathways))


results.dir = "../../results/"
snps.dir    = "../../Supplementary/" 

# load snps list
snps <- read.table(paste0(snps.dir,"SignificativeSnps_Invididual_Integrated.txt"), header = T)

# make sure all the snps are into the rsid format
if( sum(substr(snps$SNP.rsid,1,2) != "rs") ){
  stop("Not all snps are into rsid")}

# make sure the list does not cantain duplicates
snps = unique(snps$SNP.rsid)

# retrieve snps information using NCBI's refSNP
snps.to.gene = ncbi_snp_query(snps)

# add NA values in the empty cell of Gene to have a consistent data frame
snps.to.gene$gene = gsub("^$|^ $", NA, snps.to.gene$gene)


getPPI <- function(gene){
  if( is.na(gene) ) return (NA) 
  res = GET(paste0("https://string-db.org/api/tsv/network?identifiers=", gene))
  con = rawToChar(res$content)
  df  = read.table(text = con, sep = "\t", header = TRUE)
  proteins = paste( unique(df$preferredName_B), collapse = "/" )
  return (proteins)
}

# add the PPI information from STRING
snps.to.gene$PPI = sapply(snps.to.gene$gene, getPPI)

# retreive the latest GMT file listing all the genes in each of the human pathways
wp.hs.gmt = rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")

# Now we can process it to generate the two dataframes we need for enricher
wp2gene   = readPathwayGMT(wp.hs.gmt)
wpid2gene = wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name = wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

# we perform enrichment analysis from WikiPathways
wikienrich <- function(genes.entrez){
  clusterProfiler::enricher(
    genes.entrez[[2]],
    pAdjustMethod = "fdr",
    pvalueCutoff  = 0.1,   # p.adjust cutoff; relaxed because we want just the pathway name 
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
}

getPathway <- function(gene){
  if( is.na(gene) ) return (NA)
  
  # convert the SYMBOL into ENTREZID
  genes.entrez <- clusterProfiler::bitr(gene,
                                        fromType = "SYMBOL",toType = "ENTREZID",
                                        OrgDb = org.Hs.eg.db)
  
  # get the pathways set if the gene can be mapped
  if( is.null( wikienrich(genes.entrez) ) ){ 
    print(paste0("Gene: ", gene, " can't be mapped!"))
    return (NA) }
  
  ewp <- wikienrich(genes.entrez)
  
  pathways = paste0(1:nrow(ewp@result),
                    ". ",
                    ewp@result$Description, 
                    collapse = " / " )
  
  return (pathways)
}

# add the Pathways information 
snps.to.gene$Pathway = sapply(snps.to.gene$gene, getPathway)

# store the snps's information table 
path = paste0(results.dir,"SNPs_gene_information.csv")
write.table(snps.to.gene, path, quote = FALSE, row.names = FALSE, sep="\t")

