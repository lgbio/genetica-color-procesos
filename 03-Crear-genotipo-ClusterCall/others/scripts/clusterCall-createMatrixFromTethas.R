#!/usr/bin/Rscript
USAGE="
Create a matrix of Markers vs Samples from fitpoly results
INPUT: Filtered fitPoly results: MarkerName, SampleName, Theta
       start and end markers (integers), label of output 
OUTPUT: File with a Matrix of Markers vs Samples
"

source ("lglib02.R")
library (parallel)

createMatrix <- function (snp, fitThetas) {
	message ("...", snp)
	dataSnp = fitThetas [fitThetas$MarkerName==snp,]
	thetas = data.frame (t(dataSnp$Theta))
	if (all (is.na (thetas))) 
		return (NULL)
	names (thetas) = dataSnp$SampleName
	rownames (thetas) = snp
	return (thetas)
 }

#fitResultsFile = "results/out-fitPoly_scores.dat"
#fitResults = read.table (fitResultsFile, header=T, sep="\t")

#fitThetas     = fitResults [,c("MarkerName","SampleName","Theta")]
fitThetasFile = "datasets/CCC_Andigena_677_2015_fitTetra.csv"
#write.csv (fitThetas, fitThetasFile, quote=F, row.names=F)

#----
args = commandArgs(trailingOnly = TRUE)
#start = args [1]
#end   = args [2]
#label = args [3]

fitThetas  = read.csv (fitThetasFile)		
snpList    = levels (fitThetas$MarkerName)

label = "ALL"
start = 1
end   = length (snpList)

outs = mclapply (snpList[start:end], createMatrix, fitThetas, mc.cores=7)
fitThetasDF = do.call (rbind.data.frame, outs)

outFilename = addLabel ("out-ClusterCall-matrixThetas.csv", label)

fitPolyMatrix = data.frame (Makers=rownames (fitThetasDF), fitThetasDF)
write.csv (fitPolyMatrix, outFilename, quote=F, row.names=F)


