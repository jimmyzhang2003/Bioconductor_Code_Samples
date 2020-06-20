# Bioconductor Code Samples
Experimental code for genomic data science in R. Exploratory analysis of hospital data, gene set analysis of bottomly dataset, eQTL analysis of MatrixEQTL data, and DESeq2 analysis of zebrafish RNA-Seq data.

## Getting Started
### Dependencies
- R 3.6.1
- ggplot2
- Bioconductor (for package installation)

### Install Packages
Install Bioconductor: 

```r 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

Install packages from Bioconductor:

```r 
BiocManager::install(c("Biobase", "limma", "goseq", "DESeq2", "org.Mm.eg.db", "MatrixEQTL", "zebrafishRNASeq"))
```
## Usage

### Hospitals Exploratory Analysis
- download the zip file of hospital data ("https://d396qusza40orc.cloudfront.net/rprog%2Fdata%2FProgAssignment3-data.zip"), unzip the file, and read the "outcome-of-care-measures.csv" file into R
- the 11th, 17th, and 23rd columns list the mortality rates of heart attack, heart failure, and pneumonia, respectively
- the 7th column lists the state each hospital is located in

### Bottomly Data Gene Set Analysis
- extract the phenotype data, expression data, and feature data from the bottomly data 
- create a model using limma adjusting for strain and lane number, and then find p-values
- create a DESeqDataSet adjusting for strain, perform a DESeq analysis, and then calculate enrichment statistics by using ```goseq()``` on the differentially expressed results 

### Matrix eQTL 
- read in the expression, snps, and covariates files from the "MatrixEQTL" package
- adjust variables in a slice for genotype and expression data for the the file reading procedure
- match all the parameters and arguments for the Matrix eQTL Engine

### Zebrafish DESeq Analysis
- the "zebrafishRNASeq" package contains data from an zebrafish RNA-Seq experiment as a data frame called zfGenes
- the zfGenes data frame contains 92 rows with spikein transcripts, each starting with "ERCC" and should be excluded from analysis
