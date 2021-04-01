# function to create the Manhattan Plot
# R packages required to run the script:
# - tidyverse
# - ggplot2
# - ggrepel

if(!require("tidyverse", character.only = TRUE))
{
  install.packages("tidyverse")
  if(!require("tidyverse", character.only = TRUE))
  {
    stop("tidyverse package not found")
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

if(!require("ggrepel", character.only = TRUE))
{
  install.packages("ggrepel")
  if(!require("ggrepel", character.only = TRUE))
  {
    stop("ggrepel package not found")
  }
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

ManhattanGenerator <- function(gwasResults,path,pheno){
  
  df <- gwasResults %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_annotate  = ifelse(-log10(P)>=5, "yes", "no"))  %>%
    mutate( is_highlight = ifelse(is_annotate == "yes", "yes", "no")) %>%
    mutate( is_super_highlight = ifelse(-log10(P)>=7.3, "yes", "no"))  
  
  axisdf <- df %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  p = ggplot(df, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("#28286F", "#758FC5"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(limits=c(0,10),breaks=c(0,2.5,5,7.3,10))+
    # Add highlighted points
    geom_point(data=subset(df, is_highlight=="yes"), color="orange", size=2) +
    geom_point(data=subset(df, is_super_highlight=="yes"), color="red", size=3) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(df, is_annotate=="yes"), aes(label=SNP), size=2) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_text(size=13),
      axis.title.y = element_text(size=13, hjust = 0.5),
      axis.text = element_text(size=13, hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ) + ggtitle(sprintf("%s - Manhanttan plot", pheno)) + 
      xlab("Chromosome") + ylab(expression(-log[1*0](p)))
  
  p2 = p + geom_hline(yintercept=-log10(1e-5), linetype="dashed", 
                      color = "blue4", size=0.3) + 
    geom_hline(yintercept=-log10(5e-8), linetype="dashed", 
               color = "red",   size=0.3)
  ggsave(
    filename=path, 
    plot=p2, 
    device="png",
    dpi=500, 
    width = 14, 
    height = 8, 
    units = "in"
)}