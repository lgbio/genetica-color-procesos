#!/usr/bin/Rscript

#source ("lglib09.R")
suppressMessages (library (dplyr))

"Annotate selectes SNPs using PGSC genes annotations" 

#--------------------------------------------------------
# Main
#--------------------------------------------------------
main <- function () {
	message (">>> Reading inputs...")
	# INPUTS:
	#annotatedSNPs = read.csv ("inputs/annotations-SPUD-DB-Solcap-SNPs.csv", sep="\t")[,1:2]
	PGSCGenes    = read.csv ("inputs/PGSC_DM_V403_representative_genes-BGIGenes.gff", sep="\t", header=F)
	selectedSNPs = read.csv ("inputs/scores-colorsPotato-SELECTED.csv")

	# OUTPUTS:
	outFileAnnotations = "scores-colorPotato-SELECTED-ANNOTATED-SpudDB.csv"


	message (">>> Making annotations...")
	makeAnnotationsFromPGSCGenes (selectedSNPs, PGSCGenes)
}

#--------------------------------------------------------
# Make annotations
#--------------------------------------------------------
makeAnnotationsFromPGSCGenes <- function (selectedSNPs, PGSCGenes) {
	genes = PGSCGenes [,c(4,5,9)]
	colnames (genes) = c("INI", "END", "ANNOTATION")
	print (tail (genes));quit()

	funSearch <- function (POSITION) {
		gene = which (POSITION >= genes$INI & POSITION <= genes$END) 
		message ("Pos: ", POSITION, ",", genes)
		if (length (gene)) strsplit (as.character (genes [gene[1], 3]), split="name=")[[2]][1] else NA
	}

	ANN = unlist (sapply (selectedSNPs$POSITION, funSearch))
	annotations = data.frame (selectedSNPs, ANNOTATION=ANN) 

	#annotations = left_join (selectedSNPs, reprGenes)
	createDir ("outputs")
	write.csv (annotations, outFileAnnotations, quote=T, row.names=F)
}

#--------------------------------------------------------
# Make annotations
#--------------------------------------------------------
makeAnnotationsFromSpudDB <- function () {
	# Select columns and change column names
	colnames (annotatedSNPs) = c("SNP", "ANNOTATION")
	annotatedSNPs$SNP = gsub ("solcap_snp_", "", annotatedSNPs$SNP)

	annotations = left_join (selectedSNPs, annotatedSNPs)
	createDir ("outputs")
	write.csv (annotations, outFileAnnotations, quote=F, row.names=F)
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





