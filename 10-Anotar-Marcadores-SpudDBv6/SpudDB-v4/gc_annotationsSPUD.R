#!/usr/bin/Rscript

#' Annotate SNPs 
#'
#' Annotate SNPs using SPUD database annotations" 
#'
#' @param markersFile     Filename with markers (SNP column= to annotate).
#' @param annotationsFile Filename of annotations (SNP and ANNOTATION columns).
#' @return Return a dataframe with markers annotated. 
#' @export
gc_annotateMarkers <- function (markersFile, outFile, annotationsFile) {
	#makeAnnotationsFromPGSCGenes (selectedSNPs, PGSCGenes)
	makeAnnotationsFromSpudDB (markersFile, outFile, annotationsFile)
}
#--------------------------------------------------------
main <- function () {
	source ("lglib14.R")
	library (dplyr)
	# Main
	# INPUTS:
	message (">>> Reading inputs...")
	annotationsFile = "inputs/annotations-SPUD-DB-Solcap-SNPs.csv"
	PGSCGenes       = "inputs/PGSC_DM_V403_representative_genes-BGIGenes.gff"
	markersFile     = "inputs/scores-colorsPotato-SELECTED.csv"

	# OUTPUTS:
	outFile = "scores-colorPotato-SELECTED-ANNOTATED-SpudDB.csv"
	gc_annotateMarkers (markersFile, outFile, annotationsFile) 
}
#--------------------------------------------------------

#--------------------------------------------------------
# Make annotations
#--------------------------------------------------------
makeAnnotationsFromSpudDB <- function (markersFile, outFile, annotationsFile) {
	annotatedSNPs = read.csv (annotationsFile)
	view (annotatedSNPs)
	markers       = read.csv (markersFile)
	view (markers)

	message (">>> Select columns and change column names...")
	annotatedSNPs$SNP = gsub ("solcap_snp_", "", annotatedSNPs$SNP)

	annotations = left_join (markers, annotatedSNPs)

	# Move columns to first positions
	SNP        = annotations$SNP
	ANNOTATION = annotations$ANNOTATION
	annotations$SNP        = NULL
	annotations$ANNOTATION = NULL
	annotations = cbind (SNP, ANNOTATION, annotations)

	write.csv (annotations, outFile, row.names=F)
	return (annotations)
}

#--------------------------------------------------------
# Make annotations
#--------------------------------------------------------
makeAnnotationsFromPGSCGenes <- function (selectedSNPs, PGSCGenes) {
	genes = PGSCGenes [,c(4,5,9)]
	colnames (genes) = c("INI", "END", "ANNOTATION")

	funSearch <- function (POSITION) {
		gene = which (POSITION >= genes$INI & POSITION <= genes$END) 
		if (length (gene)) strsplit (as.character (genes [gene[1], 3]), split="name=")[[2]][1] else NA
	}

	ANN = unlist (sapply (selectedSNPs$POSITION, funSearch))
	annotations = data.frame (selectedSNPs, ANNOTATION=ANN) 

	#annotations = left_join (selectedSNPs, reprGenes)
	createDir ("outputs")
	write.csv (annotations, outFileAnnotations, quote=T, row.names=F)
}


#----------------------------------------------------------
# Utility for create dir, if it exists, it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", dirname (newDir), basename (newDir))
			if (dir.exists (oldDir) == T) checkOldDir (oldDir)
			file.rename (newDir, oldDir)
		}
	}
	checkOldDir (newDir)
	system (sprintf ("mkdir %s", newDir))
}
#--------------------------------------------------------
#--------------------------------------------------------
main ()


