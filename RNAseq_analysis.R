#!/usr/bin/R

args <- commandArgs(trailingOnly=TRUE)

i <- args[2]
i <- unlist(strsplit(i,split = "_"))

if(! "data.table" %in% rownames(installed.packages())){
  cat("data.table has not been installed....\nInstalling data.table\n")
  install.packages("data.table")
}

if((!"RSQLite" %in% rownames(installed.packages())) | packageDescription("RSQLite")$Version=="1.1-2"){
  cat("WARNING: Updated RSQLite version is needed... Installing updated RSQLite version...")
  install_packages("RSQLite")
}

if(! "DESeq2" %in% rownames(installed.packages())){
  cat("DESeq2 has not been installed....\nInstalling DESeq2\n")
  install.packages("DESeq2")
}


library(data.table)
library(DESeq2)

# sampleTable <-fread(args[1], header = FALSE,stringsAsFactors = TRUE)
# setnames(sampleTable,old=names(sampleTable),new=c("samplename","filename","condition"))


sampleTable <-fread(args[1], header = TRUE,stringsAsFactors = TRUE)
setnames(sampleTable,old=names(sampleTable)[1:2],new=c("samplename","filename"))


directory<-""
design_formula <- as.formula(paste("~",args[4]))


dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = design_formula) #design= ~ condition)

dds_original <- dds

# Exploratory analysis

ppca<-NULL

if("ea" %in% i ){
  dds <- dds[ rowSums(counts(dds)) > 1, ]

  rld <- rlog(dds, blind = FALSE) # blind = FALSE means that differences between cell lines and treatment (the variables in the design)
                                  # will not contribute to the expected variance-mean trend of the experiment.
  
  sampleDists <- dist(t(assay(rld)))
  
  library("pheatmap")
  library("RColorBrewer")
  
  plot_name <- gsub('DESeq_table.txt','sclust.pdf',args[1])

  ppca <- plotPCA(rld)

}

if(!is.null(ppca)){                                             # Though this if condition may seem
  pdf(file = plot_name,width=6,height=4)                        # redundant, it is necessary to save the pca plot.
}                                                               # For some weird reason, the pdf function does not
                                                                # save the pca plot inside the previous if condition.
if("ea" %in% i ){
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
               clustering_distance_rows = sampleDists,
               clustering_distance_cols = sampleDists,
               col = colors)
}
ppca

if(!is.null(ppca)){
  dev.off()
}



# Compute TPMs or FPKMs

if("tpm" %in% i) {
  library("GenomicFeatures")
  #txdb <- makeTxDbFromGFF("/home/paulafp/Documents/temp/WS255_WBcel235/c_elegans.PRJNA13758.WS255.canonical_geneset.gtf")
  txdb <- makeTxDbFromGFF(args[3])
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  
  exons.list.per.gene <- reduce(exons.list.per.gene)
  exons.list.per.gene <- as.data.table(exons.list.per.gene)
  gene_width <- exons.list.per.gene[,sum(width),by=group_name]
  
  tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
  }
  
  read_counts <- as.data.frame(counts(dds_original))
  samples <- colnames(read_counts)
  
  read_counts$group_name <- row.names(read_counts)
  read_counts <-  as.data.table(read_counts)
  setkey(read_counts,group_name)
  
  setkey(gene_width,group_name)
  
  read_counts <- read_counts[gene_width]
  
  sapply(samples,function(x){
    read_counts[,eval(x):=tpm(get(x),V1)]
  })
  
  read_counts <- as.data.frame(read_counts)
  #row.names(read_counts)<-read_counts$group_name
  read_counts$geneName <- read_counts$group_name
  read_counts$group_name <- NULL
  read_counts <- read_counts[,c("geneName",samples)]
  
  tpm_file <-gsub('DESeq_table.txt','tpm.txt',args[1])
  
  write.table(x = read_counts,file = tpm_file,quote = FALSE,sep="\t",row.names = FALSE)
  
}

# Differential Calling

if("deseq" %in% i){
  dds <- DESeq(dds)
  if(length(args)==4){
    res <- results(dds)
  }
  else{
    res <- results(dds,contrast=c(args[4:6]))
  }
  
  res <- res[order(res$padj),]
  
  print(mcols(res, use.names = TRUE))
  cat("\n")
  summary(res,alpha=0.05)
  
  res_file = gsub('DESeq_table.txt','DESeq_results.txt',args[1])
  
  res <- as.data.frame(res)
  names(res) <- c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")
  res$geneName <- rownames(res)
  res <- res[,c("geneName","baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]
  
  write.table(x = res,file = res_file,quote = FALSE,sep="\t",row.names = FALSE)
  
}


sessionInfo_file <-gsub('DESeq_table.txt','sessionInfo.txt',args[1])
writeLines(capture.output(sessionInfo()), sessionInfo_file)
