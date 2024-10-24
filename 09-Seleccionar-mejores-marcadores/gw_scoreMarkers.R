#!/usr/bin/Rscript 

#' Score GWAS markers from file
#'
#' Score and sort markers resulting from MultiGWAS by using an own score called "AgroScore". 
#' The AgroScore takes into account three terms: Genomic Control (80%), replicability (10%) and significance (10%). 
#'
#' @param scoresFile Filename of scores resulting from GWAS.
#' @param outFile    Filename of output markers with AgroScores.
#' @return None. It writes results in a file.
#' @export
gw_scoreMarkersFile <- function (scoresFile, outFile="out-ScoredMarkers.csv") {
	scores     = read.csv (scoresFile);view(scores)
	scoresAgro = gw_scoreMarkers (scores)
	write.csv (scoresAgro, outFile , row.names=F)
	return (outFile)
}

#' Score GWAS markers from dataframe
#'
#' Score and sort markers resulting from MultiGWAS by using an own score called "AgroScore". 
#' The AgroScore takes into account three terms: Genomic Control (80%), replicability (10%) and significance (10%). 
#'
#' @param scoresFile Filename of scores resulting from GWAS.
#' @param outFile    Filename of output markers with AgroScores.
#' @return dataframe with results.
#' @export
gw_scoreMarkers <- function (scores) {
	# For normalization, get number of tools, traits, componentes
	nTools  = length (unique (scores$TOOL))
	nTraits = length (unique (scores$TRAIT))
	MAXREP  = nTools * nTraits * 3

	# Replicability score: Count of SNPs between all SNPs
	scoresRepl  = scores %>% add_count (SNP, sort=F, name="ScoresRepl");#view (scoresRepl)
	valuesRepl  = scoresRepl$ScoresRepl/MAXREP;#view (valuesRepl)

	# Significance score: 1 for significants, 0 otherwise
	valuesSign   = ifelse (scores$SIGNIFICANCE, 1,0)
	scoresSign   = cbind (scoresRepl, ScoresSign=valuesSign);#view (scoresSign)

	# GC score: Measures closenes of Genomic Control (GC) to 1
	valuesGC     = 1 - abs (1-scores$GC)
	scoresGC     = cbind (scoresSign, ScoreGC=valuesGC);#view (scoresGC)

	# Difference score: Measures difference between threshold and Score 
	# It's not useful if it isn't normalized by each tool
	#valuesDiff   = scores$SCORE - scores$THRESHOLD
	#scoresDiff   = cbind (scoresGC, ScoreDiff=valuesDiff);#view (scoresDiff)

	totalScore = 0.7*valuesGC + 0.2*valuesSign + 0.1*valuesRepl;#view (valuesRepl)

	scoresGC$GSCORE = round (totalScore, 3)
	scoresAgro = scoresGC %>% arrange (desc(GSCORE))

	return (scoresAgro)
}
#-------------------------------------------------------------
# Get base trait or stem word from full trait ("PULPA.HLC.C" --> "PULPA")
#-------------------------------------------------------------
addBaseTrait <- function (scores) {
	getTrait <- function (trait) {
		strsplit (trait, split="[.]")[[1]][1]
	}
	mainTraits = sapply (as.character (scores$TRAITHCL), getTrait)
	newScores = data.frame (TRAIT=mainTraits, scores)
	return (newScores)
}

#----------------------------------------------------------- 
# Score SNPs according to GC, DIFF, REPL
#----------------------------------------------------------- 
calculateScoreSelection <- function (scoreTable) {
	nrm <- function (x, minX="", maxX="") {
		minX = if (minX=="") min(x) else minX
		maxX = if (maxX=="") max(x) else maxX
		return ( (x - minX)/(maxX-minX))
	}	

	scoreGC   = 1 - abs (1-scoreTable$GC)
	scoreGC   = sapply (scoreGC, function (x) if (x>0) x else 0.0)
	scoreREPL = nrm (scoreTable$REPL, 1, 4)
	scoreDIFF = nrm (scoreGC*scoreTable$DIFF)
	scoreSIGN = sapply (scoreTable$SIGNIFICANCE, function(x) if (x==T) 1.0 else 0.0)

	scoreSel  = 0.5*scoreGC + 0.25*scoreREPL + 0.1*scoreDIFF + 0.15*scoreSIGN

	scoreTableSNP = data.frame (scoreTable, scoreSel, scoreGC, scoreREPL, scoreDIFF, scoreSIGN)
	return (scoreTableSNP)
}

#----------------------------------------------------------- 
# Main
#----------------------------------------------------------- 
main <- function () {
	#!/usr/bin/Rscript
	suppressMessages (library (dplyr))
	source ("lglib14.R")
	options (width=300, warn=2)

	args = commandArgs(trailingOnly = TRUE)
	

	#scoresFile = "inputs/ADDITIVE-BESTMARKERS-ALLTRAITS.csv"
	scoresFile = args [1]
	message ("Scoring markers from file ", scoresFile)
	gw_scoreMarkers (scoresFile)
}

#----------------------------------------------------------- 
#----------------------------------------------------------- 
#main ()


