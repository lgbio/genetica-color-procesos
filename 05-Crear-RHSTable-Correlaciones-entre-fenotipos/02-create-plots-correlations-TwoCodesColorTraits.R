#!/usr/bin/env Rscript
USAGE="
Correlations and plots for both two coding systems and base color traits
USAGE: correlate-ColorsCodifications.R <Accessions with two codes for each phenotype>"

# Libraries
suppressMessages (library (dplyr))
library (rcompanion)
library (ggplot2)
source ("lglib14.R")
library(corrplot)

# Input/output files
INPUTFILE = "outputs/out-colors-TwoCodes-Phenos.csv"
OUTFILE01 = "outputs/out-ColorsCorrelations-TwoCodes.pdf"
OUTFILE02 = "outputs/out-ColorsCorrelations-AllPhenos.pdf"

main <- function () {
		data  = read.csv (INPUTFILE)
		data  = data [,-1]

		# Correlations between base color phenotypes
		correlationsBasePhenos (data);

		# Correlations between pairs of codifications
		correlationsCodingPhenos (data)
}

# Calculate correlations between pair of base and derived
correlationsCodingPhenos <- function (data) {
		nCols    = ncol (data)
		colorsDf = data.frame ()
		for (i in seq (1, nCols, 2)) {
			colorsTable = data [,i:(i+1)]
			traitName   = names (colorsTable)[2]
			corr        = calcCorrelation (colorsTable)
			df          = data.frame (traitName, corr)
			colorsDf    = rbind (colorsDf, df)
			print (corr)
		}
		names (colorsDf) = c ("TRAIT", "CORRELATION")
		write.csv (colorsDf, "outputs/out-ColorsCorrelations-TwoCodes.csv", row.names=F)
		CORRTEXTS = paste0 (colorsDf$CORRELATION, "%")

		print (colorsDf)

		ggplot (colorsDf, aes (x=TRAIT, y=CORRELATION, fill=TRAIT)) +
			geom_bar (stat="identity") +
			ggtitle ("Correlations between trait values for two different years") +
			geom_text(aes(label=CORRTEXTS), position=position_dodge(width=0.9), vjust=-0.25)

		ggsave (OUTFILE01)
}

# Calculate correations between multiple variables
correlationsBasePhenos <- function (twoCodesTable) {
	N         = ncol (twoCodesTable)
	# Table of RHS Codes
	RHSTable  = twoCodesTable [, seq (2,N,2)]; #view (RHSTable)
	N         = ncol (RHSTable)
	corrMatrix = matrix (nrow=N, ncol=N)
	for (x in 1:N) {
		for (y in 1:N) {
			if (is.na (corrMatrix [x,y])) {
				xydf = RHSTable [,c(x,y)]; #view (xydf) 
				corrMatrix [x,y] = calcCorrelation (xydf)
			}
		}
	}
	corrDF = data.frame (corrMatrix, row.names=names (RHSTable)); #view (corrDF)
	names (corrDF) = names (RHSTable)
	write.csv (corrDF, "outputs/out-ColorsCorrelations-AllPhenos.csv")

	# Plot the correlation matrix
	mat = as.matrix (corrDF/100); #view (mat)

	pdf (OUTFILE02)

	corrplot (mat, method = "number", type="upper", bg="black", is.corr=F, mar=c(0,0,1,0),
			  title="Correlations between potato color traits")
	dev.off()
}

# Calculate correlation between pair of variables
calcCorrelation <- function (colorsTable) {
	colnames  = names (colorsTable)
	IDVAR     = colnames [2]
	TIMEVAR   = colnames [1]
	outFile   = sprintf ("colorCorrelations-matrix-%s.csv", IDVAR)

	colorsTable = colorsTable [complete.cases (colorsTable),]
	counts = count (colorsTable, across(all_of(IDVAR)), across(all_of(TIMEVAR)))
	counts = rename (counts, Freq = n)

	# Create wide table
	countsWide = reshape (counts, idvar = IDVAR, timevar = TIMEVAR, direction="wide")
	names (countsWide) = c ("RHS", gsub ("Freq[.]", "s", names (countsWide[,-1])))
	countsWide = data.frame (countsWide, row.names=1)
	countsWide [is.na (countsWide)] = 0
	#write.csv (countsWide, outFile, row.names=F)

	# Calc categorical correlations
	countsMat = as.matrix (countsWide)
	countsMat [,-1] = as.integer (countsMat [,-1])
	corr = round (100*cramerV (countsMat), 0)
	return (corr)
}

main ()

