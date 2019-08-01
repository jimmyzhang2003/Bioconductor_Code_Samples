library(Biobase)
library(MatrixEQTL)

base.dir = find.package("MatrixEQTL")
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep = "")
expression_file_name = paste(base.dir, "/data/GE.txt", sep = "")
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep = "")
outputfile_name = tempfile()

#read in files
expr = read.table(expression_file_name, sep = "\t", header = TRUE, row.names = 1)
snps = read.table(SNP_file_name, sep = "\t", header = TRUE, row.names = 1)
cvrt = read.table(covariates_file_name, sep = "\t", header = TRUE, row.names = 1)

#parameters
pvOutputThreshold = 1e-2
errorCovariance = numeric()
useModel = modelLINEAR

#genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 1000
snps$LoadFile(SNP_file_name) 

#repeat for expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 1000
gene$LoadFile(expression_file_name)

#empty covariates object
cvrt = SlicedData$new()

#Matrix eQTL
eqtlmat = Matrix_eQTL_engine(
          snps = snps,
          gene = gene, 
          cvrt = cvrt,
          output_file_name = NULL,
          pvOutputThreshold = pvOutputThreshold,
          useModel = useModel,
          errorCovariance = errorCovariance,
          verbose = TRUE,
          pvalue.hist = TRUE,
          min.pv.by.genesnp = FALSE,
          noFDRsaveMemory = FALSE)

#analysis
plot(eqtlmat)
eqtlmat$all$eqtls

