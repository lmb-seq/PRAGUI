library(data.table)

args <- commandArgs(trailingOnly=TRUE)

tpm <- fread(args[1],header=T)

colSums(tpm[,names(tpm)[2:ncol(tpm)],with=FALSE])
