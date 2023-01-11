library(DESeq2)

countData  <- read.csv("20220530_miRNA_count_in_dme_CLASH.csv",header=T, row.names=1, sep=",")
colData <- DataFrame(condition=factor(c('ctrl','ctrl','ctrl','z8','z8','z8')))
ddsHTSeq <- DESeqDataSetFromMatrix(countData, colData,design = ~ condition)
dds <- DESeq(ddsHTSeq, betaPrior = FALSE) # differential expression
res <- results(dds, contrast=c("condition","z8","ctrl"))
write.csv(res,file='20220530CLASH_ctrl_z8.csv', sep= ", ", row.names=TRUE,col.names=TRUE)

