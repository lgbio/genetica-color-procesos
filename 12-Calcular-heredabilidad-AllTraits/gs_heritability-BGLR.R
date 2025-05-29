#!/usr/bin/env Rscript

#------------------- For calling from a renv environment -------------------------------------------------
Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa/14x-SeleccionGenomica-N-Marcadores")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
#---------------------------------------------------------------------------------------------------------

NCORES = 8
DEBUG  = FALSE
if (DEBUG==FALSE) {
	INPUTSDIR="inputs/"           # Real training dir
}else { # Only for testing
	INPUTSDIR="inputs-test/"     # Toy datasets, only for testing the script
}
#-------------------------------------------------------------
# Main only for testing
#-------------------------------------------------------------
main <-  function() {
	#library (ppGS)
	source ("lglib14.R")
	source ("gs_plotsHCL.R")
	library (parallel)
	library (dplyr)
	library (ggplot2)
	library (BGLR)

    phenotypesFile   = paste0 (INPUTSDIR, "fenotiposColor-ComponentesLCH.csv")
    genotypeFile     = paste0 (INPUTSDIR, "genotipo-AndigenaCCC-ClusterCall2020-MATRIX.csv")
    outputDir        = "outputs"

    hclTableFilename = gs_heritability (genotypeFile, phenotypesFile, outputDir)

    plotHeritabilities (hclTableFilename)
}

#'-------------------------------------------------------------

#' Calculates heritabilities for HCL penotypes.
#' 
#' @param genoFile   Filename of genotype matrix.
#' @param phenosFile Filename of phenotypes table.
#' @param outputDir  Dirname where outputs will be saved.
#' @return None. Generate two files: a table (CSV) with heritabilities for each trait
#'         and a summary plot (PDF) comparing heritabilities from all traits.
#' @import parallel
#' @import dplyr
#' @import ggplot2
#' @import BGLR
#' @export
gs_heritability <- function (genoFile, phenosFile, outputDir) {
	data       = readProcessGenoPheno (genoFile, phenosFile)
	geno       = data$geno
	pheno      = data$pheno 
	traitNames = colnames (pheno)

	message (">>> Calculating heritabilities...")
	createDir (outputDir)
	createDir (sprintf ("%s/tmp", outputDir))

	h2Results = mclapply (traitNames, calculateHeritabilityTrait, 
						  geno, pheno, outputDir, mc.cores=NCORES)
	hclTableFilename   = createHCLTable (traitNames, h2Results, outputDir)

    return (hclTableFilename)
}

#-------------------------------------------------------------
# Calculata heritability for a single trait using BayesA method
#-------------------------------------------------------------
calculateHeritabilityTrait <- function (traitName, geno, pheno, outputDir) {
	message (">>> Heritability for ", traitName)
	# Get Trait and Geno for BGLR
	X = scale(geno)/sqrt(ncol(geno))
	Y = pheno [, traitName]
	 
	# Estimation using variance components
	MODEL  = "BayesA"
	PREFIX = paste0 (outputDir, "/tmp/", MODEL, "-", traitName,"-")
	#varUFilename = paste0 (PREFIX, "ETA_1_varB.dat")
	varUFilename = paste0 (PREFIX, "ETA_1_ScaleBayesA.dat")
	varEFilename = paste0 (PREFIX, "varE.dat")

	eta  = list(list(X=X,model=MODEL,saveEffects=T))
	fm   = BGLR(y=Y,ETA=eta, nIter=15000, burnIn=5000, thin=5, verbose=F, saveAt=PREFIX)
	varU = scan(varUFilename)
	varE = scan(varEFilename)
	h2   = varU/(varU+varE)

	return (h2)
}

#-------------------------------------------------------------
# Plot heritabilities dot, bar, and boxes
#-------------------------------------------------------------
plotHeritabilities <- function (hclTableFilename) {
    # Bar Plot (Recommended for Clarity)

    # Assuming your data is in a data frame called 'data'
    # with columns "Traits" and "Heritability"

	data  = read.csv (hclTableFilename, check.names=F); #view (data)

    # Create the bar plot to PDF
    ggplot (data, aes(x = reorder(Traits, -Heritability), y = Heritability)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      coord_flip() +  # Flips axes for better readability
      labs(title = "Heritability of Traits", x = "Traits", y = "Heritability") +
      theme_minimal()

	outPlotname   = gsub ("TABLE.csv", "PLOT-BAR.pdf", hclTableFilename)
	ggsave (outPlotname, width=11)

    # Saving Dot Plot to PDF
    dot_plot <- ggplot(data, aes(x = reorder(Traits, Heritability), y = Heritability)) +
      geom_point(size = 4, color = "blue") +
      coord_flip() +
      labs(title = "Heritability of Traits", x = "Traits", y = "Heritability") +
      theme_minimal()

	outPlotname   = gsub ("TABLE.csv", "PLOT-DOT.pdf", hclTableFilename)
	ggsave (outPlotname, width=11)

    # Saving boxplots to PDF
	gs_plotsHCL (hclTableFilename, "outputs", "Heritability", 
				 "Comparison of mean heritabilities for phenotypes")
}

#-------------------------------------------------------------
# Plot heritability by multi boxplots
#-------------------------------------------------------------
createHCLTable <- function (traitNames, h2Results, outputDir) {
	message (">>> Creating HCL table...")

	h2Table = data.frame ()
	for (i in 1:length (h2Results)) {
		Prefix        = strsplit (traitNames[i], "[.]")[[1]][1]
		HCL_component = strsplit (traitNames[i], "[.]")[[1]][2]
		Models        = "BA"
		Traits        = traitNames [i]
		Heritability  = h2Results [[i]]
		tmpDF  = data.frame (Prefix, HCL_component, Models, 
							 Traits, Heritability,stringsAsFactors=F)
		h2Table = rbind (h2Table, tmpDF)
	}
	hclTable = h2Table [order (h2Table$Traits),]
	hclTableFilename  = sprintf ("%s/out-%s-HCL-Comparison-TABLE.csv", outputDir, "Heritability")
	write.csv (hclTable, hclTableFilename , quote=F, row.names=F)
	return (hclTableFilename)
}


#-------------------------------------------------------------
# Read genotype and phenotype, transform, impute, and filter by MAF
#-------------------------------------------------------------
readProcessGenoPheno <- function (genotypeFile, phenotypeFile) {
	# Read data, get geno matrix, and process it for BGRL
	genoCCC  = read.csv (genotypeFile, check.names=F, row.names=1);#view(genoCCC)
	phenoCCC = read.csv (phenotypeFile, check.names=F,  row.names=1);#view (phenoCCC)
 
	genoCCC  = t (genoCCC); #view (genoCCC)
   
	samplesGeno   = rownames (genoCCC);#view (samplesGeno)
	samplesPheno  = rownames (phenoCCC);#view (samplesPheno)
	samplesCommon = intersect (samplesGeno, samplesPheno); #view (samplesCommon)
	view (samplesCommon)

	geno  = genoCCC  [samplesCommon,]; #view (geno)
	pheno = phenoCCC [samplesCommon,]; #view (pheno)

	view (geno)
	view (pheno)

	#-------------------------------------------------------------
	# Impute NA alleles
	#-------------------------------------------------------------
	imputeGenotype <- function (M) {
		message(">>> Missing marker data imputed with population mode...")
		impute.mode <- function(x) {
			ix <- which(is.na(x))
			if (length(ix)>0) 
				x[ix] <- as.integer(names(which.max(table(x))))
			return(x)
		}
		missing <- which(is.na(M))
		if (length(missing)>0) {
			M <- apply(M,2,impute.mode)
		}
		return (M)
	}
	genoImputed = imputeGenotype (geno);#view (genoImputed)

	#-------------------------------------------------------------}
	# MAF genotype
	#-------------------------------------------------------------}
	MAFGenotype <- function (M, thresholdMAF) {
		message (">>> Checking minor allele frecuency, MAF=", thresholdMAF)
		calcMAF <- function(x) {
			ploidy=4
			AF <- mean(x, na.rm=T) / ploidy;
			MAF <- ifelse(AF > 0.5,1-AF,AF)
		}
		MAF <- apply(M, 2, calcMAF)
		polymorphic <- which(MAF>thresholdMAF)
		M <- M[,polymorphic]
	}
	genoMAF = MAFGenotype (genoImputed, 0.1);#view(genoMAF)

	return (list(geno=genoMAF, pheno=pheno))
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
	dir.create (sprintf (newDir))
}

#-------------------------------------------------------------
# Call main
#-------------------------------------------------------------
main ()

