#!/usr/bin/R

args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
setwd(wd)

if(! "data.table" %in% rownames(installed.packages())){
  cat("data.table has not been installed....\nInstalling data.table\n")
  install.packages("data.table")
}

if(! "devtools" %in% rownames(installed.packages())){
  cat("devtools has not been installed....\nInstalling devtools\n")
  install.packages("devtools")
}
library(devtools)

if((!"RSQLite" %in% rownames(installed.packages())) | packageDescription("RSQLite")$Version!="1.1-2"){
  cat("WARNING: RSQLite version 1.1-2 is needed for cummeRbund... Installing RSQLite version 1.1-2...")
  install_version("RSQLite", version = "1.1-2", repos = "http://cran.us.r-project.org")
}

if(!"cummeRbund" %in% rownames(installed.packages())){
  cat("cummeRbund has not been installed....\nInstalling cummeRbund\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("cummeRbund")
}


library(data.table)
library(cummeRbund)

cuff<-readCufflinks()

s<-csScatterMatrix(genes(cuff))

dend<-csDendro(genes(cuff),replicates = T)

genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)

gene.diff<-diffData(genes(cuff))

gene.features<-annotation(genes(cuff))

gene_diff2 <- gene.diff[abs(gene.diff$log2_fold_change)>1 & gene.diff$significant=="yes" &
                          abs(gene.diff$log2_fold_change) != Inf ,]
gene_diff2$abs_log2_fold_change <- abs(gene_diff2$log2_fold_change)
gene_diff2 <-gene_diff2[order(gene_diff2$abs_log2_fold_change,decreasing = TRUE),]
myGeneIds <- gene_diff2$gene_id[0:50]


myGeneIds <- unique(gene.features[gene.features$gene_id %in% myGeneIds,]$gene_short_name)
myGenes<-getGenes(cuff,myGeneIds)

h<-csHeatmap(myGenes,cluster='both')


pdf(file = "exploratory_analysis_plots.pdf")
# Scattermatrix
s
# Dendrogram
plot(dend)
# MDS plot
genes.MDS.rep
# Heatmap
h
dev.off()

