#!/usr/bin/env Rscript

#------------------- For calling from a renv environment -------------------------------------------------
Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa/14x-SeleccionGenomica-N-Marcadores")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
#---------------------------------------------------------------------------------------------------------


#------------ Running parameters -----------------------------
NCORES = 8
DEBUG  = FALSE
if (DEBUG == FALSE) {
	INPUTSDIR="inputs/"			  # Real training dir
}else {  # Only for testing
	INPUTSDIR="inputs-test/"	 # Toy datasets, only for testing the script
}

#----------- Model Parameters --------------------------------
NITER    = 15000
BURNIN   = 5000
THIN     = 5
VERBOSE  = FALSE
MODEL    = "BRR"      # BRR, BayesA, RKHS

ADDVARFILESUFFIX = switch (MODEL,
		BayesA = "ScaleBayesA.dat",
		BRR    = "varB.dat",
		RKHS   = "varU.dat"
		)

#-------------------------------------------------------------
# Main only for testing
#-------------------------------------------------------------
main <-  function() {
	#library (ppGS)
	source ("lglib14.R")
	source ("gs_plotsHCL.R")

	library (BGLR)
	library (parallel)
	library (ggplot2)
	library (dplyr)
	library(AGHmatrix)

	phenotypesFile	 = paste0 (INPUTSDIR, "fenotiposColor-ComponentesLCH.csv")
	genotypeFile	 = paste0 (INPUTSDIR, "genotipo-AndigenaCCC-ClusterCall2020-MATRIX.csv")
	outputDir		 = paste0 ("outputs-", MODEL)

	hclTableFilename = gs_heritability (genotypeFile, phenotypesFile, outputDir)
    plotHeritabilities (hclTableFilename, outputDir)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
gs_heritability <- function (genoFile, phenosFile, outputDir) {
	data	   = readProcessGenoPheno (genoFile, phenosFile)
	geno	   = data$geno
	view (geno)
	pheno	   = data$pheno 
	traitNames = colnames (pheno)

	message (">>> Calculating heritabilities...")
	#createDir (outputDir)
	#createDir (sprintf ("%s/tmp", outputDir))

	GmatrixPoly = calculate_Gmatrix (geno)
	view (GmatrixPoly)

	h2Results = mclapply (traitNames, heritabilityPolyTrait, 
						  geno, pheno, GmatrixPoly,  outputDir, mc.cores=NCORES)

	hclTableFilename   = createHCLTable (traitNames, h2Results, outputDir)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
heritabilityPolyTrait <- function (traitName, geno, pheno, GmatrixPoly, outputDir) {
	# Output dirs and filenames
	OUTFILEPREFIX   = paste0 (outputDir, "/tmp/", MODEL, "-", traitName,"-")
	additiveVarFile = paste0 (OUTFILEPREFIX, "ETA_1_", ADDVARFILESUFFIX)    # Additive variance for RKHS
	residualVarFile = paste0 (OUTFILEPREFIX, "varE.dat")          # Residual variance for RKHS

	if (file.exists (additiveVarFile) && file.exists (residualVarFile)) {
		var_additive = scan(additiveVarFile)
		var_residual = scan(residualVarFile)
		h2_narrow    = var_additive / (var_additive + var_residual)
		return (h2_narrow)
	}
		
	# Extract the phenotype column (e.g., 'CTallo.L') and assign IDs from 'Registro'
	phenoVector = pheno [, traitName]
	if (MODEL == "RKHS") {
		genoMatrix  = GmatrixPoly
		ETA         = list (list (K = genoMatrix, model = MODEL))
	}else {  # BayesA, BRR
		genoMatrix = scale (geno)/sqrt(ncol(geno))
		ETA        = list (list(X = genoMatrix, model = MODEL))
	}

	fit <- BGLR( y = phenoVector, ETA = ETA,
	  nIter = NITER, burnIn = BURNIN, thin = THIN , verbose = VERBOSE, saveAt = OUTFILEPREFIX 
	)

	var_additive = scan(additiveVarFile)
	var_residual = scan(residualVarFile)
	h2_narrow    = var_additive / (var_additive + var_residual)

	return (h2_narrow)
}

#--------------------------------------------------------------------------------
# Assuming your genotype data is `genotype_matrix_0_4` (individuals x markers)
# coded as 0, 1, 2, 3, 4 for a tetraploid.
#--------------------------------------------------------------------------------
calculate_Gmatrix <- function (geno) {
	message ("+++ Calculating Gmatrix...")
	genotype_matrix_0_4 = geno
	ploidy_level <- 4 # For tetraploid
	GmatrixPoly <- AGHmatrix::Gmatrix(
	  SNPmatrix = genotype_matrix_0_4,
	  method = "VanRaden", # A common method; "Slater" is another option for polyploids
	  #method = "Slater", # A common method; "Slater" is another option for polyploids
	  ploidy = ploidy_level,
	  ploidy.correction = T
	  # Other arguments like missingValue, maf, etc., as needed
	)
	return (GmatrixPoly)
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
		Models        = MODEL
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
# Plot heritabilities dot, bar, and boxes
#-------------------------------------------------------------
plotHeritabilities <- function (hclTableFilename, outputDir) {
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
	gs_plotsHCL (hclTableFilename, outputDir, "Heritability", 
				 "Comparison of mean heritabilities for phenotypes")
}

#-------------------------------------------------------------
# Read genotype and phenotype, transform, impute, and filter by MAF
#-------------------------------------------------------------
readProcessGenoPheno <- function (genotypeFile, phenotypeFile) {
	# Read data, get geno matrix, and process it for BGRL
	genoCCC  = read.csv (genotypeFile, check.names=F, row.names=1);#view(genoCCC)
	phenoCCC = read.csv (phenotypeFile, check.names=F,	row.names=1);#view (phenoCCC)
 
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
		
#--------------------------------------------------------------------
#--------------------------------------------------------------------
main ()
