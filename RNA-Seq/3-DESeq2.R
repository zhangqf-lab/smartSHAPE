
library("DESeq2")

infn0 <- '/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/5.cufflinks/genecount_0h.csv'
infn3 <- '/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/5.cufflinks/genecount_3h.csv'

genecount0 <- read.csv(infn0, header=TRUE); head(genecount0)
genecount3 <- read.csv(infn3, header=TRUE); head(genecount3)

condition <- c('WT','WT','WT','KO','KO','KO')
coldata <- data.frame(condition, row.names=c('WT_1','WT_2','WT_3','KO_1','KO_2','KO_3'))

rownames(genecount0) <- as.list(genecount0['X0'])[[1]]
genecount0 <- genecount0[ c('WT_1','WT_2','WT_3','KO_1','KO_2','KO_3') ]

rownames(genecount3) <- as.list(genecount3['X0'])[[1]]
genecount3 <- genecount3[ c('WT_1','WT_2','WT_3','KO_1','KO_2','KO_3') ]

dds0 <- DESeqDataSetFromMatrix(countData=genecount0,colData=coldata, design=~condition)
dds3 <- DESeqDataSetFromMatrix(countData=genecount3,colData=coldata, design=~condition)

keep0 <- rowSums(counts(dds0)) >= 10
dds0 <- dds0[keep0,]
dds0 <- DESeq(dds0)
res0 <- results(dds0)

keep3 <- rowSums(counts(dds3)) >= 10
dds3 <- dds3[keep3,]
dds3 <- DESeq(dds3)
res3 <- results(dds3)

resOrdered0 <- res0[order(res0$pvalue),]; summary(res0)
resOrdered3 <- res3[order(res3$pvalue),]; summary(res3)

plotMA(res0, ylim=c(-2,2), main='DESeq2 0h')
plotMA(res3, ylim=c(-2,2), main='DESeq2 3h')

write.csv(resOrdered0, 
    file="/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/5.cufflinks/DESeq2_0h.csv",
    quote=FALSE)
write.csv(resOrdered3, 
    file="/Share2/home/zhangqf7/meiling/data/low_input_icshape/2019-12-24-ko/Processing/5.cufflinks/DESeq2_3h.csv",
    quote=FALSE)

# Cd274: ENSMUSG00000016496.7
plotCounts(dds0, gene='ENSMUSG00000016496.7', intgroup="condition", main="Cd274 0h")
plotCounts(dds3, gene='ENSMUSG00000016496.7', intgroup="condition", main="Cd274 3h")

# Il6: ENSMUSG00000025746.11
plotCounts(dds0, gene='ENSMUSG00000025746.11', intgroup="condition", main="Il6 0h")
plotCounts(dds3, gene='ENSMUSG00000025746.11', intgroup="condition", main="Il6 3h")

# Ptgs2: ENSMUSG00000032487.8
plotCounts(dds0, gene='ENSMUSG00000032487.8', intgroup="condition", main="Ptgs2 0h")
plotCounts(dds3, gene='ENSMUSG00000032487.8', intgroup="condition", main="Ptgs2 3h")