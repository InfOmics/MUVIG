library(magrittr)
library(clusterProfiler)

#  ****************************** download the WikiPathways ******************************
wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) 
wpid2name <- wp2gene %>% dplyr::select(wpid, name) 

setwd("../Without_Swedd/")
files = setdiff(list.files(), list.dirs(recursive = FALSE, full.names = FALSE))

for (i in files){
name = substring(i,10)
name = gsub("\\..*","", name)

table <- read.csv(paste(i,sep=""))
################## PREPARE INPUT ##################### 
geneList = data.frame(table$ENTREZID,2^table$logFC)
geneList = na.omit(geneList)
geneList = geneList[!duplicated(geneList$table.ENTREZID),]
geneList <- sort(deframe(geneList), decreasing = T)
genes = names(geneList)[abs(geneList) > 1]
################## CORE PART ##################### 
ewp  <- enricher(genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name) 
ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE) 
# convert the gene IDs to gene symbols
ewp  <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
ewp2 <- setReadable(ewp2, org.Hs.eg.db, keyType = "ENTREZID")

write.csv(ewp@result,paste("Enricher/",i,"_enricher.csv",sep=""), quote=F,row.names = F)
write.csv(ewp2@result,paste("GSEA/",i,"_enricher.csv",sep=""), quote=F,row.names = F)
}
