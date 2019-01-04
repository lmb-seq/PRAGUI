library(data.table)
library(DESeq2)
library(refGenome) # Provides functions for reading and summarising GTF files

organism = "human"

databases <- list(human=list("org.Hs.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","REFSEQ","SYMBOL","UCSCKG")),
               mouse=list("org.Mm.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","MGI","REFSEQ","SYMBOL")),
               worm=list("org.Ce.eg.db",c("ACCNUM","ALIAS","ENSEMBL","ENTREZID","GENENAME","MGI","REFSEQ","SYMBOL","WORMBASE")),
               fly=list("org.Dm.eg.db"), # NEEDS TO BE COMPLETED!!
               fish=list("org.Dm.eg.db"), # NEEDS TO BE COMPLETED!!
               yeast=list("org.Sc.sgd.db")) # NEEDS TO BE COMPLETED!!


if(organism %in% names(databases)){
  dtb <- databases[[organism]][[1]]
  ids <- databases[[organism]][[2]]
  names(dtb)<-NULL
  library(Biobase)
  library(AnnotationDbi)
  library(dtb,character.only = TRUE)
}

sampleTable <-fread("/data/Eszter/fastq/GSE48213/samples_cancerlines.csv_DESeq_table.txt")

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                   directory = "",
                                   design = ~ condition)
dds<-DESeq(dds)

baseMeanPerLvl <- as.data.frame(sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) ))
baseMeanPerLvl$gene_id <- rownames(baseMeanPerLvl)
baseMeanPerLvl <- as.data.table(baseMeanPerLvl)

#  Create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome() 

wd <- getwd()
setwd("/data/genome_assemblies/Homo_sapiens/")
# read GTF file into ensemblGenome object
read.gtf(ens, "Homo_sapiens.GRCh38.91.gtf")
setwd(wd)

# create table of genes
my_gene <- as.data.table(getGenePositions(ens))
my_gene[,locus:=paste0(seqid,":",start,"_",end)]
if(any(!is.na(my_gene$gene_name))){
  my_gene<- my_gene[,c("gene_id","gene_name","locus"),with=FALSE]
} else {
  my_gene<- my_gene[,c("gene_id", "locus"),with=FALSE]
  if(organism %in% names(databases)){
    test<-my_gene$gene_id[1:5]
    findannot <- FALSE
    i=1
    while(findannot==FALSE){
      annot<-ids[i]
      findannot<-any(test %in% keys(org.Hs.eg.db,keytype = annot))
      i=i+1
    }
  }
  all_names <- function(x){paste(x,sep = "_")}
  geneSymbols <- unlist(mapIds(get(dtb), keys=my_gene$gene_id, column="SYMBOL", keytype=annot,multiVals=all_names))
  geneSymbols2<-data.table(gene=geneSymbols,gene_id=names(geneSymbols))
  setkey(geneSymbols2,gene_id)
  setkey(baseMeanPerLvl,gene_id)
  baseMeanPerLvl[geneSymbols2]
}

setkey(baseMeanPerLvl,gene_id)
setkey(my_gene,gene_id)
baseMeanPerLvl<- baseMeanPerLvl[my_gene]

comparisons <- combn(levels(dds$condition),2)

i=0
for(c in 1:dim(comparisons)[2]){
  print(i)
  comp <- comparisons[,c]
  res <- results(dds,contrast = c("condition", comp))
  baseMean2Lvls <- baseMeanPerLvl[,c("gene_id","gene_name","locus",comp),with=FALSE]
  res <- res[order(res$padj),]
  
  res <- as.data.frame(res)
  names(res) <- c("baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")
  res$gene_id <- rownames(res)
  res <- res[,c("gene_id","baseMean", "log2FoldChange","lfcSE","stat","pvalue","padj")]
  
  rownames(res)<-NULL
  head(res)
  res_4_pete <- as.data.table(res[,c("gene_id", "log2FoldChange","stat","pvalue","padj")])
  
  setnames(res_4_pete,names(res_4_pete),c("gene_id","log2(fold_change)","test_stat","p_value","q_value"))
  setkey(res_4_pete,gene_id)
  
  res_4_pete<-res_4_pete[baseMean2Lvls]
  res_4_pete[,significant:=ifelse(q_value<=0.05,"yes","no")]
  res_4_pete[,test_id:=gene_id]
  res_4_pete[,sample_1:=comp[1]]
  res_4_pete[,sample_2:=comp[2]]
  res_4_pete[,status:=ifelse(!is.na(p_value),"OK","NOTEST")]
  setnames(res_4_pete,c("gene_name",comp),c("gene","value_1","value_2"))
  res_4_pete <- res_4_pete[,c("test_id","gene_id","gene","locus","sample_1", "sample_2", "status", "value_1","value_2", "log2(fold_change)", "test_stat", "p_value", "q_value", "significant"), with=FALSE]
  res_4_pete<-res_4_pete[order(p_value)]
  if(i==0){
    res_4_pete_comb<-res_4_pete
  } else {
    res_4_pete_comb<-rbind(res_4_pete_comb,res_4_pete)
  }
  i<-i+1
}

write.table(x = res_4_pete_comb,file ="~/Desktop/ex_output_4_peat.txt",quote = FALSE,sep="\t",row.names = FALSE)

