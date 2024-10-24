#!/usr/bin/Rscript
source ("lglib14.R")

info="
Create HCL component phenotypes from potato phenotypes
associated with color and measured using the RHS scale
INPUTS:
	1. Table with RHS codes and HCL components
	2. Individuals with potato phenotypes associated with color

OUTPUTS:
	1. Matrix of HCL component phenotypes for each potato phenotype
"
options (width=300, stringsAsFactors=F)

#--- Inputs / Outputs -------------------------------------
inHCLTable       = "inputs/RoyalHorticultural-ColourCharts.csv"
inRHSValues      = "inputs/out-colors-RHSValues-BaseTraits.csv"
inBaseTraitNames = "inputs/out-colors-BaseShortNames.csv"
inBaseTraitNames = "inputs/trait-color-ColumNames.csv"

outPotatoHCLPhenos = "outputs/out-fenotiposColor-ComponentesLCH.csv"
outSummaryTraits   = "outputs/out-summary-TraitNames.csv"
createDir ("outputs")

#----------------------------------------------------------
# Main
#----------------------------------------------------------
main <- function () {
	source ("lglib14.R")
	# Read input data
	hclComponents   = read.csv (inHCLTable, check.names=F, row.names=1)
	RHSTraitValues  = read.csv (inRHSValues, check.names=F, header=T)
	baseTraitNames  = read.csv (inBaseTraitNames, check.names=F)

	# Get trait full and short names
	baseTraitsFullNames  = baseTraitNames [,1]
	baseTraitsShortNames = baseTraitNames [,2]

	# Create table with short names
	traitsPotato    = RHSTraitValues

	hclNames   = colnames (hclComponents)[-1:-7];#view (hclNames)
	traitNames = colnames (traitsPotato)[-1];#view (traitNames)

	hclTraitsTable = data.frame (traitsPotato[,1]) 
	colnamesList = c("Registro")

	# Summary table showing base trait and derived trait names
	summaryValues = c()
	
	for (i in 1:length(traitNames)) {
		trait = traitNames [i]
		print (trait)
		fullName  = baseTraitsFullNames [i]
		summaryValues = c (summaryValues, fullName, trait)
		for (hcl in hclNames) {
			rhsCodes       = traitsPotato [,c(trait)] 
			hData          = hclComponents [rhsCodes, hcl] 
			colName        = paste0 (trait,".",gsub ("LCH_","",hcl))
			colnamesList   = c (colnamesList, colName)
			hclTraitsTable = data.frame (hclTraitsTable, hData);#view (hclTraitsTable)

			summaryValues = c (summaryValues, colName)
		}
	}
	# Create summary table
	summaryTable = data.frame (matrix (summaryValues, ncol = 5, byrow=T ))
	names (summaryTable) = c("Base trait name", "Short name", "HCL luminance name", "HCL chroma name", "HCL hue name")

	# Write outputs
	write.csv (summaryTable, outSummaryTraits, row.names=F)
	colnames (hclTraitsTable) = colnamesList
	write.csv (hclTraitsTable, outPotatoHCLPhenos, row.names=F, quote=F)
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
#----------------------------------------------------------
main ()
