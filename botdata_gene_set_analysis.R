library(BioBase)
library(limma)
library(goseq)
library(DESeq2)
library(org.Mm.eg.db)

con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata = pData(bot)
fdata = featureData(bot)
edata = exprs(bot)
fdata = fdata[rowMeans(edata) > 5]     ###remove lowly expressed genes
edata = edata[rowMeans(edata) > 5,]

#analysis with limma
mod = model.matrix(~pdata$strain + pdata$lane.number)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)
limma_pvals = topTable(ebayes_limma, number = dim(edata)[1])$P.Value
hist(limma_pvals, col = 3)

#gene set analysis with goseq
expr = edata
grp = pdata$strain
pdf = data.frame(grp)
de = DESeqDataSetFromMatrix(expr, pdf, ~ grp)
de_fit = DESeq(de)
de_results = results(de_fit)
genes = as.integer(de_results$padj < 0.05)
not_na = !is.na(genes)
names(genes) = rownames(expr)
genes = genes[not_na]
pwf = nullp(genes, "mm9", "ensGene")
GO.wall = goseq(pwf, "mm9", "ensGene")
head(Go.wall)

