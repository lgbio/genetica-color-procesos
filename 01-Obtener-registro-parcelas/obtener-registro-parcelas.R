#!/usr/bin/Rscript

USAGE="
Get Parcela and Registro for Andigena dataParcelaRegistro"

options (stringsAsFactors=F)

fileGenotipar = "original/2019_campo_xgenotipar_BGV_421.csv"
outFilename   = "outputs/tabla-ParcelaRegistro-andigena.csv"

#-------------------------------------------------------------
# Main
#-------------------------------------------------------------
main <- function () {
	dataParcelaRegistro = read.csv (fileGenotipar)[,c(1,2)]
	colnames (dataParcelaRegistro) = c("Parcela", "Registro")

	parcelas     = dataParcelaRegistro$Parcela
	parcelasAnd  = as.character (parcelas [grepl ("And", parcelas)])
	dataAnd      = dataParcelaRegistro [dataParcelaRegistro$Parcela %in% parcelasAnd,]
	rownames (dataAnd) = dataAnd$Parcela
	dataNoNum    = dataAnd [which(is.na (as.numeric(dataAnd$Registro))),]
	dataNum      = dataAnd [setdiff (dataAnd$Parcela, dataNoNum$Parcela),]
	dataDups     = dataNum [which(duplicated (dataNum$Registro)),] 


	createDir ("outputs")
	dataErrors   = rbind (dataNoNum, dataDups)
	write.csv (dataErrors, addLabel (outFilename, "ERRORS"), quote=F, row.names=F)

	dataCleaned  = dataAnd [setdiff (dataAnd$Parcela, dataErrors$Parcela),]
	write.csv (dataCleaned, addLabel (outFilename, "CLEANED"), quote=F, row.names=F)
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



