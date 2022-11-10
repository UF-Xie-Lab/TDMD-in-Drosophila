source("https://bioconductor.org/biocLite.R")
library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)

setwd('/orange/mingyi.xie/luli/Exp198_drosophila_CLASH_3rd')
#countData  <- read.csv("20220530_miRNA_count_in_dme_CLASH.csv",header=T, row.names=1, sep=",")
countData  <- read.csv("test.csv",header=T, row.names=1, sep=",")
#colData <- DataFrame(condition=factor(c('ctrl','ctrl','ctrl','z8','z8','z8')))
colData <- DataFrame(condition=factor(c(z8,'ctrl',z8,'z8',ctrl,ctrl)))
ddsHTSeq <- DESeqDataSetFromMatrix(countData, colData,design = ~ condition)
ddsHTSeq
#write.csv(assay(ddsHTSeq),file='a.csv', sep= ", ", row.names=TRUE,col.names=TRUE)##output array
head(assay(ddsHTSeq))#Take a look of ddsHTSeq
dds <- DESeq(ddsHTSeq, betaPrior = FALSE) # differential expression
#res <- results(dds, contrast=c("condition","ctrl","zs8"))
res <- results(dds)
#write.csv(res,file='20220530CLASH_ctrl_z8.csv', sep= ", ", row.names=TRUE,col.names=TRUE)
write.csv(res,file='test_out.csv', sep= ", ", row.names=TRUE,col.names=TRUE)
