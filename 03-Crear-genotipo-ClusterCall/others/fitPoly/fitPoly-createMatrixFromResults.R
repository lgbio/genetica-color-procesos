#!/usr/bin/Rscript
USAGE="
Create a matrix of Markers vs Samples from fitpoly results
INPUT: Filtered fitPoly results: MarkerName, SampleName, geno
       start and end markers (integers), label of output 
OUTPUT: File with a Matrix of Markers vs Samples
"

source ("lglib02.R")
library (parallel)

createMatrix <- function (snp, fitGenos) {
	message ("...", snp)
	dataSnp = fitGenos [fitGenos$MarkerName==snp,]
	genos = data.frame (t(dataSnp$geno))
	if (all (is.na (genos))) 
		return (NULL)
	names (genos) = dataSnp$SampleName
	rownames (genos) = snp
	return (genos)
 }

fitGenosFile = "outputs/out-predictions-CCC-fitPoly-TABLE.csv"

#----
args = commandArgs(trailingOnly = TRUE)

fitGenos  = read.csv (fitGenosFile)		
snpList    = levels (fitGenos$MarkerName)

outs = mclapply (snpList, createMatrix, fitGenos, mc.cores=7)
fitGenosDF = do.call (rbind.data.frame, outs)

outFilename = "outputs/out-predictions-CCC-fitPoly-MATRIX.csv"

fitPolyMatrix = data.frame (Makers=rownames (fitGenosDF), fitGenosDF)
write.csv (fitPolyMatrix, outFilename, quote=F, row.names=F)

# Invert genotypes (0..4, 1..3, 2..3, 3..1, 4..0)
inFile         = outFilename
dataMatrix     = read.csv (inFile)
invertedMatrix = data.frame (Markers=dataMatrix [,1], 4-dataMatrix[,-1])

newFile = gsub (".csv", "-INVERTED.csv", inFile)
write.csv (invertedMatrix, newFile, quote=F, row.names=F)

