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
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2")
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

if(!"tximport" %in% packages){
  cat("tximport has not been installed....\nInstalling data.table\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("tximport")
}

if(!"GenomicFeatures" %in% packages){
  cat("GenomicFeatures has not been installed....\nInstalling data.table\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",ask=FALSE)
  BiocManager::install("GenomicFeatures")
}

databases <- list(human=list("org.Hs.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","REFSEQ","SYMBOL","UCSCKG")),
                  mouse=list("org.Mm.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","MGI","REFSEQ","SYMBOL")),
                  worm=list("org.Ce.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","REFSEQ","SYMBOL","WORMBASE")),
                  fly=list("org.Dm.eg.db", c("ACCNUM","ALIAS","ENSEMBL", "ENTREZID","FLYBASE","GENENAME","REFSEQ","SYMBOL")), 
                  zebrafish=list("org.Dr.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","REFSEQ","SYMBOL","ZFIN")),
                  yeast=list("org.Sc.sgd.db",c("ALIAS","ENSEMBL","ENTREZID","GENENAME","REFSEQ","SGD"))) 

gtf <- args[3]
organism <- args[4]



library(data.table)
library(DESeq2)
library(GenomicFeatures)



# sampleTable <-fread(args[1], header = FALSE,stringsAsFactors = TRUE)
# setnames(sampleTable,old=names(sampleTable),new=c("samplename","filename","condition"))


sampleTable <-fread(args[1], header = TRUE,stringsAsFactors = TRUE)
setnames(sampleTable,old=names(sampleTable)[1:2],new=c("samplename","filename"))


directory<-""
design_formula <- as.formula(paste("~",args[5]))

txdb <- makeTxDbFromGFF(gtf)

if("salmon" %in% i){
  library(tximport)
  files <- sampleTable$filename
  files <- as.character(files)
  names(files) <- sampleTable$samplename
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- select(txdb, k, "GENEID", "TXNAME")
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  colnames(txi$counts)<-sampleTable$samplename
  sampleTable <- as.data.frame(sampleTable)
  sampleTable2<-data.frame(sampleTable[,args[5]])
  colnames(sampleTable2)<-args[5]
  rownames(sampleTable2)<-sampleTable$samplename
  dds <- DESeqDataSetFromTximport(txi, sampleTable2, design_formula)
  } else {
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = design_formula) #design= ~ condition)
  }
  
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
  
  ppca <- plotPCA(rld,intgroup=args[5])

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
  tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
  }
  
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  
  exons.list.per.gene <- reduce(exons.list.per.gene)
  exons.list.per.gene <- as.data.table(exons.list.per.gene)
  gene_width <- exons.list.per.gene[,sum(width),by=group_name]
  
  read_counts <- as.data.frame(counts(dds_original))
  samples <- colnames(read_counts)

  read_counts$group_name <- row.names(read_counts)
  read_counts <-  as.data.table(read_counts)
  setkey(read_counts,group_name)

  setkey(gene_width,group_name)

  read_counts <- read_counts[gene_width]
  
  read_counts <- read_counts[!is.na(get(names(read_counts)[1]))]
  
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

print(dim(combn(levels(dds[[args[5]]]),2)))

if("deseq" %in% i){
  dds <- DESeq(dds)
  
  nc <- as.data.frame(counts(dds,normalized=TRUE))
  nc$gene_id <- rownames(nc)
  nc<- nc[,c(colnames(nc)[length(colnames(nc))],colnames(nc)[-length(colnames(nc))])]
  
  nc_file <-gsub('_table.txt','_norm_read_counts.txt',args[1])
  
  write.table(x = nc,file = nc_file,quote = FALSE,sep="\t",row.names = FALSE)
  
  baseMeanPerLvl <- as.data.frame(sapply( levels(dds[[args[5]]]), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds[[args[5]]] == lvl] ) ))
  baseMeanPerLvl$gene_id <- rownames(baseMeanPerLvl)
  baseMeanPerLvl <- as.data.table(baseMeanPerLvl)
  print(head(baseMeanPerLvl))
  
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
  my_gene[,locus:=paste0(seqid,":",start,"_",end)]
  if(any(!is.na(my_gene$gene_name))){
    my_gene<- my_gene[,c("gene_id","gene_name","locus"),with=FALSE]
  } else {
    my_gene<- my_gene[,c("gene_id", "locus"),with=FALSE]
    if(organism %in% names(databases)){
      dtb <- databases[[organism]][[1]]
      ids <- databases[[organism]][[2]]
      names(dtb)<-NULL
      inst_pkgs <- c()
      if(!"Biobase" %in% packages){
        cat("Biobase has not been installed...\n")
        inst_pkgs <- c(inst_pkgs,"Biobase")
        }
      if(!"AnnotationDbi" %in% packages){
        cat("AnnotationDbi has not been installed...\n")
        inst_pkgs <- c(inst_pkgs,"Biobase")
        }
      if(!dtb %in% packages){
        cat(paste0(dtb," has not been installed...\n"))
        inst_pkgs <- c(inst_pkgs,dtb)
        }
      if(length(inst_pkgs)>0){
        cat(paste0("Installing ",inst_pkgs, "...\n"))
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(inst_pkgs,ask = FALSE)
        }
      library(Biobase)
      library(AnnotationDbi)
      library(dtb,character.only = TRUE)
      test <- my_gene$gene_id
      findannot <- FALSE
      i=1
      while(findannot==FALSE){
        annot<-ids[i]
        findannot<-any(test %in% keys(get(dtb),keytype = annot))
        i=i+1
      }
      all_names <- function(x){paste(x,sep = "_")}
      geneSymbols <- unlist(mapIds(get(dtb), keys=my_gene$gene_id, column="SYMBOL", keytype=annot,multiVals=all_names))
      geneSymbols2<-data.table(gene_name=geneSymbols,gene_id=names(geneSymbols))
      setkey(geneSymbols2,gene_id)
      setkey(baseMeanPerLvl,gene_id)
      baseMeanPerLvl <- baseMeanPerLvl[geneSymbols2]
    } else {
      baseMeanPerLvl[,gene_name:=NA]
    }
  }
  
  setkey(baseMeanPerLvl,gene_id)
  setkey(my_gene,gene_id)
  baseMeanPerLvl<- baseMeanPerLvl[my_gene]
  
  #comparisons <- combn(levels(dds[[args[5]]]),2)
  
  if(length(args)==5){
    # comparisons <- combn(levels(dds[[args[5]]]),2)
    comparisons <- combn(levels(dds[[args[5]]]),2)
  }
  else{
    comparisons <- matrix(args[6:7],nrow = 2)
  }
  
  i=0
  for(c in 1:dim(comparisons)[2]){
    print(i)
    comp <- comparisons[,c]
    res <- results(dds,contrast = c(args[5], comp))
    baseMean2Lvls <- baseMeanPerLvl[,c("gene_id","gene_name","locus",comp),with=FALSE]
    res <- res[order(res$padj),]
    
    res <- as.data.frame(res)
    names(res) <- c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")
    res$gene_id <- rownames(res)
    res <- res[,c("gene_id","baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]
    
    rownames(res)<-NULL
    head(res)
    res_4_peat <- as.data.table(res[,c("gene_id", "log2FoldChange","stat","pvalue","padj")])
    
    setnames(res_4_peat,names(res_4_peat),c("gene_id","log2(fold_change)","test_stat","p_value","q_value"))
    setkey(res_4_peat,gene_id)
    
    res_4_peat<-res_4_peat[baseMean2Lvls]
    res_4_peat[,significant:=ifelse(q_value<=0.05,"yes","no")]
    res_4_peat[,test_id:=gene_id]
    res_4_peat[,sample_1:=comp[1]]
    res_4_peat[,sample_2:=comp[2]]
    res_4_peat[,status:=ifelse(!is.na(p_value),"OK","NOTEST")]
    setnames(res_4_peat,c("gene_name",comp),c("gene","value_1","value_2"))
    res_4_peat <- res_4_peat[,c("test_id","gene_id","gene","locus","sample_1", "sample_2", "status", "value_1","value_2", "log2(fold_change)", "test_stat", "p_value", "q_value", "significant"), with=FALSE]
    res_4_peat<-res_4_peat[order(p_value)]
    if(i==0){
      res_4_peat_comb<-res_4_peat
    } else {
      res_4_peat_comb<-rbind(res_4_peat_comb,res_4_peat)
    }
    i<-i+1
  }
  res_file <-gsub('table.txt','results_4_peat.txt',args[1])
  
  write.table(x = res_4_peat_comb,file = res_file,quote = FALSE,sep="\t",row.names = FALSE)
  
}

sessionInfo_file <-gsub('DESeq_table.txt','sessionInfo.txt',args[1])
writeLines(capture.output(sessionInfo()), sessionInfo_file)
