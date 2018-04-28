#!/usr/bin/R

args <- commandArgs(trailingOnly=TRUE)

i <- args[2]
i <- unlist(strsplit(i,split = "_"))

lmb_clust_packages = "/lmb/home/paulafp/applications/" # Only needed to run pipeline at the LMB

if("R_lib" %in% dir(lmb_clust_packages)){
  pack_loc = paste0(lmb_clust_packages,"R_lib/")
  packages = rownames(installed.packages(pack_loc))
  .libPaths(c(.libPaths(),pack_loc))
} else {
  packages = rownames(installed.packages())
}

if(!"data.table" %in% packages){
  cat("data.table has not been installed....\nInstalling data.table\n")
  install.packages("data.table",repos='http://cran.us.r-project.org')
}

if((!"RSQLite" %in% packages) | packageDescription("RSQLite")$Version=="1.1-2"){
  cat("WARNING: Updated RSQLite version is needed... Installing updated RSQLite version...")
  install.packages("RSQLite",repos='http://cran.us.r-project.org')
}

if(!"DESeq2" %in% packages){
  cat("DESeq2 has not been installed... Installing DESeq2\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
}

if(!"pheatmap" %in% packages){
  cat("pheatmap has not been installed....\nInstalling data.table\n")
  install.packages("pheatmap",repos='http://cran.us.r-project.org')
}

if(!"RColorBrewer" %in% packages){
  cat("RColorBrewer has not been installed....\nInstalling data.table\n")
  install.packages("RColorBrewer",repos='http://cran.us.r-project.org')
}

if(!"refGenome" %in% packages){
  cat("refGenome has not been installed....\nInstalling data.table\n")
  install.packages("refGenome",repos='http://cran.us.r-project.org')
}


databases <- list(human=list("org.Hs.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","REFSEQ","SYMBOL","UCSCKG")),
                  mouse=list("org.Mm.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","MGI","REFSEQ","SYMBOL")),
                  worm=list("org.Ce.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","MGI","REFSEQ","SYMBOL","WORMBASE")),
                  fly=list("org.Dm.eg.db"), # NEEDS TO BE COMPLETED!!
                  zebrafish=list("org.Dm.eg.db"), # NEEDS TO BE COMPLETED!!
                  yeast=list("org.Sc.sgd.db")) # NEEDS TO BE COMPLETED!!

gtf <- args[3]
organism <- args[4]


if(organism %in% names(databases)){
  dtb <- databases[[organism]][[1]]
  ids <- databases[[organism]][[2]]
  names(dtb)<-NULL
  source("https://bioconductor.org/biocLite.R")
  biocLite(pkgs = c("Biobase","AnnotationDbi",dtb),ask = FALSE)
  library(Biobase)
  library(AnnotationDbi)
  library(dtb,character.only = TRUE)
}

library(data.table)
library(DESeq2)

# sampleTable <-fread(args[1], header = FALSE,stringsAsFactors = TRUE)
# setnames(sampleTable,old=names(sampleTable),new=c("samplename","filename","condition"))


sampleTable <-fread(args[1], header = TRUE,stringsAsFactors = TRUE)
setnames(sampleTable,old=names(sampleTable)[1:2],new=c("samplename","filename"))


directory<-""
design_formula <- as.formula(paste("~",args[5]))


dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = design_formula) #design= ~ condition)

dds_original <- dds

# Exploratory analysis

ppca<-NULL

if("ea" %in% i ){
  dds <- dds[ rowSums(counts(dds)) > 1, ]

  rld <- rlog(dds, blind = FALSE) # blind = FALSE means that the experimental design is used in estimating
                                  # the global amount of variability in the counts.
                                  # However, it is not used directly in the transformation of read counts.
  
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
  txdb <- makeTxDbFromGFF(gtf)
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
  if(length(args)==5){
    res <- results(dds)
  }
  else{
    res <- results(dds,contrast=c(args[5:7]))
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

if(organism %in% names(databases)){
  baseMeanPerLvl <- as.data.frame(sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) ))
  colnames(baseMeanPerLvl)<-paste0("sample_",1:length(colnames(baseMeanPerLvl)))
  baseMeanPerLvl$gene_id <- rownames(baseMeanPerLvl)
  baseMeanPerLvl <- as.data.table(baseMeanPerLvl)
  
  res_4_pete <- as.data.table(res[,c("geneName", "log2FoldChange","stat","pvalue","padj")])
  
  setnames(res_4_pete,names(res_4_pete),c("gene_id","log2(fold_change)","test_stat","p_value","q_value"))
  
  #  Create ensemblGenome object for storing Ensembl genomic annotation data
  library(refGenome)
  ens <- ensemblGenome() 
  
  wd <- getwd()
  setwd(dirname(gtf))
  # read GTF file into ensemblGenome object
  read.gtf(ens, basename(gtf))
  setwd(wd)
  
  # create table of genes
  my_gene <- as.data.table(getGenePositions(ens))
  my_gene[,locus:=paste0(seqid,":",start,"-",end)]
  
  setkey(res_4_pete,gene_id)
  setkey(baseMeanPerLvl,gene_id)
  setkey(my_gene,gene_id)
  
  cols <- c("gene_id","locus",names(baseMeanPerLvl)[-length(baseMeanPerLvl)],"log2(fold_change)","test_stat","p_value","q_value")
  res_4_pete<-res_4_pete[baseMeanPerLvl]
  res_4_pete<-res_4_pete[my_gene]
  res_4_pete<-res_4_pete[,cols,with=FALSE]
  
  res_4_pete[,significant:=ifelse(q_value<=0.05,"yes","no")]
  
  test<-res_4_pete$gene_id[1:5]
  findannot <- FALSE
  i=1
  while(findannot==FALSE){
    annot<-ids[i]
    findannot<-any(test %in% keys(get(dtb),keytype = annot))
    i=i+1
  }
  
  all_names <- function(x){paste(x,sep = "_")}
  
  geneSymbols <- unlist(mapIds(get(dtb), keys=res_4_pete$gene_id, column="SYMBOL", keytype=annot,multiVals=all_names))
  
  geneSymbols2<-data.table(gene=geneSymbols,gene_id=names(geneSymbols))
  setkey(geneSymbols2,gene_id)
  
  res_4_pete<-geneSymbols2[res_4_pete]
  
  res_4_pete[,test_id:=gene_id]
  
  res_4_pete<-res_4_pete[order(p_value)]
  res_4_pete <- res_4_pete[,c("test_id","gene_id","gene","locus","sample_1", "sample_2", "log2(fold_change)", "test_stat", "p_value", "q_value", "significant"), with=FALSE]
  
  res_4_pete_file = gsub('DESeq_table.txt','DESeq_results_4_pete.txt',args[1])
  write.table(x = res_4_pete,file = res_4_pete_file,quote = FALSE,sep="\t",row.names = FALSE)
}

sessionInfo_file <-gsub('DESeq_table.txt','sessionInfo.txt',args[1])
writeLines(capture.output(sessionInfo()), sessionInfo_file)
