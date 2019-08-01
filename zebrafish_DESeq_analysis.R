library(DESeq2)
library(zebrafishRNASeq)

data(zfGenes)
rownames = rownames(zfGenes)
removedrows = zfGenes[grep('^ERCC', rownames),]   ##remove spikein transcripts
zfGenes = zfGenes[!removedrows,]

#DESeq2 Analysis
countmat = as.matrix(zfGenes)
expdesign = data.frame(row.names = colnames(countmat),
                       condition = c(rep("untreated", 3), rep("treated", 3)),
                       type = c(rep("ANY", 6)))
dds = DESeqDataSetFromMatrix(countData = countmat, 
                             colData = expdesign,
                             design = ~ condition)
dds = DESeq(dds)

#results
res = results(dds)
res = na.omit(res)
res = res[order(res$padj)]
nrow(res[res$padj <= 0.05,])   ##number of features differentially expressed between control and treatment
