#!/usr/bin/Rscript
USAGE="Get genotype matrices and apply different filters:
- Select only 'andigena' accessions (And_XXX' accesions)
- Replace accesion Ids: from parcelas to register numbers
- Remove monomorfic markers (e.g 4 4 4 4 4 ..)
"
library (parallel)
#source ("lglib09.R")

main <-function () {
	registerTableFilename = "inputs/tabla-ParcelaRegistro-andigena-CLEANED.csv" 
	ratiosTableFilename   = "inputs/ratios-thetas-papa-andigena.csv"
	outputDir             = "outputs/"

	changeParcelaIdToRegisterNumber (ratiosTableFilename, registerTableFilename, outputDir)
}

#----------------------------------------------------------
# Change accessions Ids from parcelas Ids to register numbers
#----------------------------------------------------------
changeParcelaIdToRegisterNumber <- function (tableFilename, registerFilename, outputDir) {
	message ("Changing parcela to register for: ", tableFilename)
	dataWithRegs  = read.csv (registerFilename,row.names=1,check.names=F)
	dataGenos = read.csv (tableFilename)

	accessionNames = colnames (dataGenos[,-1])

	#message ("Getting registers...")
	accessionsWithRegs = c ()
	accessionsWithoutRegs = c()
	accessionsRegs = c()
	for (accession in accessionNames) {
		reg = dataWithRegs [accession,] 
		if (is.na (reg))
			accessionsWithoutRegs = c (accessionsWithoutRegs, accession)
		else {
			accessionsWithRegs = c (accessionsWithRegs, accession)
			accessionsRegs = c (accessionsRegs, reg)
		}
	}

	#message ("Changing accession names...")
	accessions = dataGenos [, accessionsWithRegs]
	colnames (accessions) = accessionsRegs

	dataWithRegs = data.frame (Markers=dataGenos [,1], accessions, check.names=F)

	createDir (outputDir)
	#message ("Writing regs table...")
	regFilename = paste0 (outputDir, addLabel (basename (tableFilename), "REGISTROS"))
	write.csv (dataWithRegs, regFilename, quote=F, row.names=F)

	colNames = c (colnames (dataGenos)[1], accessionsWithoutRegs)
	dataWithoutRegs = dataGenos [, colNames]
	#message ("Writing no regs table...")
	noregFilename = paste0 (outputDir, addLabel (basename (tableFilename), "SINREGISTROS"))
	write.csv (dataWithoutRegs, noregFilename, quote=F, row.names=F)

	return (regFilename)
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
#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}

#----------------------------------------------------------
main ()
