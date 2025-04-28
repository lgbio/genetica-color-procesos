#!/usr/bin/Rscript

INFO="
	1. It uses the ClusterCall package to generate genotype calls for autotetraploid samples from single
	   nucleotide polymorphism (SNP) marker array data. From theta values for F1 populations (such as biparental
	   families), genotype cluster positions are calibrated and with marker genotypes (training dataset)
	   for two or more F1 populations, genotypes can be predicted for a population of samples of arbitrary 
	   composition, such as a diversity panel.
	2. Remove monomorfic markers (e.g 4 4 4 4 4 ..)
	3. Add markers info to matrix
"

# Train and prediction of genotypes from theta ratios
# It uses three families to train 
# INPUTS: - Trainin datasets with six files with theta and radio files for the three families
#         - Input theta ratios for the population to predict the genotypes
# OUTPUT: - One file with predictions only for the population
# AUTHOR: Luis Garreta (2020) (lgarreta@gmail.com)

# Load libraries and set global parameters
library (parallel)
library (devtools)
load_all("lib_ClusterCall")

#source ("lglib09.R")

NCORES = 4
DEBUG  = FALSE

#---- INPUTS ----
#---- Input trainig datasets ----------
markersInfoFile = "inputs/mappinginfo-markers-Berdugo2018.csv" # Input markers information (marker|chrom|pos|ref|alt) 
#---- Input training datasets and ratios file ----
if (DEBUG==FALSE) {
	TRAINING_DIR="training-datasets/"           # Real training dir
	ratioThetasFile = "inputs/ratios-thetas-papa-andigena-REGISTROS.csv" # Real markers dataset (markers X individuals)
}else { # Only for testing
	TRAINING_DIR="training-datasets/samples50/"     # Toy datasets, only for testing the script
	ratioThetasFile = paste0 (TRAINING_DIR, "ratios-thetas-papa-andigena-REGISTROS.csv")  # Toy markers data, only for testing the script
}

training_thetaFiles = paste0 (TRAINING_DIR, c ("AxS_theta.csv", "RGxP_theta.csv", "WxL_theta.csv"))
training_radioFiles = paste0 (TRAINING_DIR, c ("AxS_r.csv", "RGxP_r.csv", "WxL_r.csv"))

#---- OUTPUT FILE ----
OUTPUTFILEGENO = "outputs/genotipo-AndigenaCCC-ClusterCall.csv"
OUTPUTFILEMATRIX = "outputs/out-predictions-clusterCall2020-MATRIX.csv"

#----------------------------------------------------------
#----------------------------------------------------------
main <- function () {
	message (">>>> Processing families....")
	trainPredictGenotypes (ratioThetasFile, training_thetaFiles, training_radioFiles, OUTPUTFILEMATRIX)
	polyFilename = removeMonomorficMarkers (OUTPUTFILEMATRIX)
	genoFilename = createGWASpolyGenotype (polyFilename, markersInfoFile)
	removePrefixSolcapSnp (genoFilename, OUTPUTFILEGENO)
}

#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
removePrefixSolcapSnp <- function (genoFilename, OUTPUTFILEGENO) {
	geno = read.csv (genoFilename, check.names=F)
	print (geno[1:5,1:5])
	geno$Markers = gsub ("solcap_snp_","", geno$Markers)
	write.csv (geno, OUTPUTFILEGENO, row.names=F, quote=F)
}

#--------------------------------------------------------------------------------------------------------
# Train and predict genotypes using ClusterCall
#--------------------------------------------------------------------------------------------------------
trainPredictGenotypes <- function (ratioThetasFile, training_thetaFiles, training_radioFiles, OUTPUTFILEMATRIX) {
	fun_readPop <- function (inputFilesItem) {
		tFile = inputFilesItem [1]
		rFile = inputFilesItem [2]
        message (">>> tFile: ", tFile)
        message (">>> rFile: ", rFile)
        theta = read.csv (tFile)
        r     = read.csv (rFile)

        cnt = colnames(theta); cnr = colnames(r); rnt = rownames(theta); rnr =  rownames(r);

        message (">>> cnt: ", length (cnt), " crn: ", length (cnr), " rnt: ", length (rnt), " rnr: ", length (rnr))

		out <- read.pop(theta.file = tFile, r.file = rFile, error.checking=F)
	}

	# Create pairs of arguments for mclapply function
	inputFiles = Map (c, training_thetaFiles, training_radioFiles)
	pops = mclapply (inputFiles, fun_readPop, mc.cores=NCORES)

	as <- pops [[1]]
	rr <- pops [[2]]
	wl <- pops [[3]]

	# Call genotypes for F1 populations
	# The CC.bipop function makes genotype calls in F1 populations:
	AS <- CC.bipop (as, parent1 = "Atlantic", parent2 = "Superior", n.core=NCORES)
	RR <- CC.bipop (rr, parent1 = "RioGrandeRusset", parent2 = "PremierRusset", n.core=NCORES)
	WL <- CC.bipop (wl, parent1 = "Wauseon", parent2 = "Lenape", n.core=NCORES)

	# Genotype an arbitrary population
	diversity <- read.pop(theta.file = ratioThetasFile, error.checking = T, thresh = 2, cex = 0.5)

	# Use the F1 populations as training data to make genotype calls for the 
	# diversity panel with the function CC.anypop.
	message (">>>> Running prediction...")
	train <- list(AS, RR, WL)
	prediction <- CC.anypop(train = train, predict = diversity[[1]], impute = T, n.core=NCORES)

	message (">>>> Writing results...")
	message (">>>>     Creating output directory...")

	createDir ("outputs")

	# Write results only for the predicted population
	colsPrediction = ncol (prediction@geno)
	colsFamilies   = ncol(as@theta) + ncol (rr@theta) + ncol (wl@theta)
	colsPopulation = (colsFamilies+1):colsPrediction # Counts two names in each subdataset

	predictionPop = prediction
	predictionPop@geno = predictionPop@geno [,colsPopulation]
	predictionPop@geno = predictionPop@geno [, order (colnames (predictionPop@geno))]

	write.pop(predictionPop, file = OUTPUTFILEMATRIX, thresh = 0.95)
}

#----------------------------------------------------------
# Filter from table monomorfic markers (e.g. 4 4 NA ... 4)
#----------------------------------------------------------
removeMonomorficMarkers <- function (regFilename) {
	message (">>>> Removing monomorfics for: ", regFilename)
	#---- local function---------------------------------------------
	checkMonomorficMarker <- function (genotypeList) {
		noNAs = Filter(function (x) !is.na(x), unlist (genotypeList))
		if (all(noNAs[1]==noNAs)) 
			return (NULL)
		return (genotypeList)
	}#--------------------------------------------------------------
	dataTable  = read.csv (regFilename, row.names=1, check.names=F)

	outMarkers = apply (dataTable, 1, checkMonomorficMarker)
	if (is.list (outMarkers))
		newMatrix = do.call (rbind, outMarkers)
	else
		newMatrix = t(outMarkers)

	# Polymorfic
	newDataTable = data.frame (Markers=rownames(newMatrix), newMatrix, check.names=F )
	polyFilename  = addLabel (regFilename,"POLYMORFIC")
	write.csv (newDataTable, polyFilename, quote=F, row.names=F)

	# Monomorfic
	monomorfic  = dataTable [setdiff (rownames (dataTable), newDataTable$Markers), ]
	monoFilename = addLabel (regFilename, "MONOMORFIC")
	write.csv (monomorfic, monoFilename, quote=F, row.names=T)

	return (polyFilename)
}

#-----------------------------------------------------------------------
# From matrix of numeric genotypes creates a GWASpoly ACGT genotype with
# chromosome and position info
#-----------------------------------------------------------------------
createGWASpolyGenotype <- function (matrixFilename, mapsFilename) {
	genos = read.csv (matrixFilename, check.names=F)
	maps  = read.csv (mapsFilename)
	genomap = merge (maps, genos, by.x="Markers", by.y="Markers", all.x=F, all.y=T)

	#---------- local fun ------------
	numbersToACGTs <- function (genomap) {
		snp     = genomap [1]
		chr     = genomap [2]
		pos     = genomap [3]
		ref     = genomap [4]
		alt     = genomap [5]
		numbers = genomap [-1:-5]

		numToACGT <- function (numberAlt, alt, ref) {
			if (is.na (numberAlt)) 
				return (NA)
			numberRef = 4 - as.numeric (numberAlt)
			return (paste0(strrep (alt, numberAlt),strrep (ref, numberRef)))
		}

		acgts = lapply (numbers, numToACGT,  alt, ref)
		rowACGTs = c (snp, chr, pos, acgts)
		return (rowACGTs)
	}
	#---------- local fun ------------
	acgtsList    = apply (genomap,1, numbersToACGTs);# acgtsList 
	acgtsMat     = do.call (rbind, acgtsList)
	genoFilename = addLabel (matrixFilename, "MAP")
	write.csv (acgtsMat, genoFilename, quote=F, row.names=F)
	return (genoFilename)
}

#----------------------------------------------------------
# Utility for create dir, if it exists, it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) {
				checkOldDir (oldDir)
			}
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
#----------------------------------------------------------
main ()

