#!/usr/bin/env Rscript

#------------------- For calling from a renv environment -------------------------------------------------
Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa/14x-SeleccionGenomica-N-Marcadores")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
#---------------------------------------------------------------------------------------------------------

#------------ Running parameters -----------------------------
NCORES = 1
DEBUG  = TRUE
if (DEBUG==FALSE) {
	INPUTSDIR="inputs/"			  # Real training dir
}else {  # Only for testing
	INPUTSDIR="inputs-test/"	 # Toy datasets, only for testing the script
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
	library(rrBLUP) # We need rrBLUP for the A.mat function
	library(BWGS)

    phenotypesFile   = paste0 (INPUTSDIR, "fenotiposColor-ComponentesLCH.csv")
    genotypeFile     = paste0 (INPUTSDIR, "genotipo-AndigenaCCC-ClusterCall2020-MATRIX.csv")
	outputDir		 = "outputs"

	hclTableFilename = gs_heritability (genotypeFile, phenotypesFile, outputDir)

	plotHeritabilities (hclTableFilename)
}

#'-------------------------------------------------------------

#' Calculates heritabilities for HCL penotypes.
#' 
#' @param genoFile	 Filename of genotype matrix.
#' @param phenosFile Filename of phenotypes table.
#' @param outputDir  Dirname where outputs will be saved.
#' @return None. Generate two files: a table (CSV) with heritabilities for each trait
#'		   and a summary plot (PDF) comparing heritabilities from all traits.
#' @import parallel
#' @import dplyr
#' @import ggplot2
#' @import BGLR
#' @export
gs_heritability <- function (genoFile, phenosFile, outputDir) {
	data	   = readProcessGenoPheno (genoFile, phenosFile)
	geno	   = data$geno
	pheno	   = data$pheno 
	traitNames = colnames (pheno)

	message (">>> Calculating heritabilities...")
	createDir (outputDir)
	createDir (sprintf ("%s/tmp", outputDir))

	#h2Results = mclapply (traitNames, calculateHeritabilityTrait, 
	#					  geno, pheno, outputDir, mc.cores=NCORES)

	#h2Results = calculateHeritabilityTrait ("CTallo.L", geno, pheno, outputDir)
	h2Results = heritabilityTrait ("CTallo.L", geno, pheno, outputDir)

	quit ()
	hclTableFilename   = createHCLTable (traitNames, h2Results, outputDir)

	return (hclTableFilename)
}

#-------------------------------------------------------------
# Calculata heritability for a single trait using BayesA method
#-------------------------------------------------------------
new_calculateHeritabilityTrait <- function(traitName, geno, pheno, outputDir) {
	print (paste0 (">>> Calculating Heritability for ", traitName))

##	# Preprocess geno and pthen
#	geno = as.data.frame (geno, check.names=F)
#	print (rownames (geno))
#	print (rownames (pheno))
##
#	pheno <- cbind (Registro = as.character (rownames(pheno)), pheno)
#	pheno = pheno [, c("Registro", traitName)]
#	view (pheno)
#	geno  <- cbind(Registro = rownames(geno),  geno)	
#	view (geno)
##
#	stopifnot (identical(rownames(geno), rownames(pheno)))
#	stopifnot (nrow(geno) == nrow(pheno))	
##
#	# Assuming you have ID columns called "ID" in both
#	merged_data <- merge(pheno, geno, by = "Registro")
#	# Then separate back into pheno and geno
#	pheno_columns = colnames (pheno)
#	geno_columns = colnames (pheno)
#	pheno <- merged_data[,pheno_columns]
#	geno <- merged_data[,geno_columns]	

	
# 1. Set 'Registro' as row names and drop the column
#rownames(pheno) <- pheno$Registro  
pheno <- pheno[, "CTallo.L", drop = FALSE]  # Keep only the phenotype column  
view (pheno)

# Extract the phenotype column (e.g., 'CTallo.L') and assign IDs from 'Registro'
pheno_vector <- pheno$CTallo.L             # Extract trait values
names(pheno_vector) <- rownames (pheno)

#rownames(geno) <- geno$Registro  
#geno <- geno[, -1]  # Drop the 'Registro' column (now in row names)  

# 2. Convert pheno to a numeric vector (required by bwgs.cv)
#pheno <- as.numeric(pheno$CTallo.L)  
#names(pheno) <- rownames(pheno)  # Preserve IDs  

# 3. Convert geno to a numeric matrix (required by bwgs.cv)
geno <- as.matrix(geno) 
view (geno)

# 4. Verify alignment
stopifnot(identical(rownames(geno), names(pheno_vector)))

		# Run bwgs.cv with 5-fold CV, 3 repetitions
		results <- bwgs.cv(
		  geno = geno,
		  pheno = pheno_vector,
		  predict.method = "RR",     # RR-BLUP method
		  geno.impute.method = "NULL",  # or "beagle" if available
		  MAF = 0.05,                # Minor allele frequency filter
		  MAXNA = 0.5,               # Remove SNPs with more than 50% missing
		  r2 = 0.0,                  # LD pruning threshold (0 means no pruning)
		  nFolds = 2,  #5,                # Cross-validation folds
		  nTimes = 1   #3                 # Number of repetitions
		)

		narrow_h2 = results$h2

		# Print narrow-sense heritability
		print (results$h2)


	# Return heritability estimate
	return (narrow_h2)
}

heritabilityTrait  <- function (traitName, geno, pheno, outputDir) {
  nIter = 15000
  burnIn = 5000
  thin = 5
  verbose = FALSE
  ploidy  = 4

  genotype = geno
  phenotype <- pheno$CTallo.L             # Extract trait values
  names(phenotype) <- rownames (pheno)
  view (phenotype)

  # 1. Input Validation
  if (!is.numeric(phenotype)) {
    stop("Phenotype must be a numeric vector.")
  }
  if (!is.matrix(genotype) || !is.numeric(genotype)) {
    stop("Genotype must be a numeric matrix.")
  }
  if (!is.numeric(ploidy) || ploidy < 1 || ploidy %% 1 != 0) {
    stop("Ploidy must be a positive integer.")
  }

  # Check for individual IDs/names - CRUCIAL for alignment
  if (is.null(rownames(genotype))) {
    stop("Genotype matrix must have row names (individual IDs) for proper alignment with phenotype.")
  }
  if (is.null(names(phenotype))) {
    stop("Phenotype vector must have names (individual IDs) for proper alignment with genotype.")
  }

  # **CRITICAL CHECK: Validate Genotype Coding based on PLOIDY**
  # Genotype values should be between 0 and `ploidy` for dosage data.
  unique_genotype_values <- unique(as.vector(genotype[!is.na(genotype)]))
  unique_genotype_values <- unique_genotype_values[!is.na(unique_genotype_values)] # Remove NAs

  if (length(unique_genotype_values) > 0) {
    if (!all(unique_genotype_values %in% 0:ploidy)) {
      # Identify problematic values
      problematic_values <- sort(unique_genotype_values[!unique_genotype_values %in% 0:ploidy])

      # Stop if any value exceeds the ploidy level (e.g., 5 for a tetraploid)
      if (any(problematic_values > ploidy)) {
         stop(paste0("Genotype data contains values greater than the specified ploidy (", ploidy, "). Problematic values: ",
                     paste(problematic_values[problematic_values > ploidy], collapse = ", "),
                     ". For allele dosage data, values should be between 0 and ploidy. Please recode your genotype matrix."))
      }
      # Warn if values are negative (though rare)
      if (any(problematic_values < 0)) {
         warning(paste0("Genotype data contains negative values: ", paste(problematic_values[problematic_values < 0], collapse = ", "), ". Please check your genotype coding."))
      }
    }
  } else {
    stop("Genotype matrix is empty or contains only NA values after filtering.")
  }


  # 2. Calculate the Genomic Relationship Matrix (G)
  message(paste0("Calculating Genomic Relationship Matrix (G) using rrBLUP::A.mat with ploidy = ", ploidy, "... This might take a moment."))
  # Pass the ploidy argument to A.mat
  G <- rrBLUP::A.mat(genotype, ploidy = ploidy)

  # Additional checks on G matrix quality
  if (any(is.na(G)) || any(is.infinite(G))) {
      stop("The genomic relationship matrix (G) generated by rrBLUP::A.mat contains NA or Inf values. This often indicates severe issues with the input genotype data (e.g., non-standard coding, too few markers, or extreme allele frequencies, or a large number of missing values).")
  }

  # Ensure G is symmetric and positive semi-definite (numerical stability check)
  if (!isSymmetric(G)) { # Check if G is symmetric
      warning("G matrix is not perfectly symmetric; forcing symmetry for numerical stability.")
      G <- (G + t(G)) / 2
  }
  if (!is.positive.semidefinite(G)) {
    epsilon <- 1e-5 # A small value for the ridge
    message("Genomic Relationship Matrix (G) is not positive semi-definite. Adding a small ridge (", epsilon, ") to the diagonal for numerical stability.")
    G <- G + diag(epsilon, nrow(G))
    if (!is.positive.semidefinite(G)) { # Check again after adding ridge
        stop("Genomic Relationship Matrix (G) remains non-positive semi-definite even after adding a ridge. This is a severe numerical stability issue; please check your genotype data for extreme collinearity or lack of variation.")
    }
  }


  # 3. Align Phenotype and G Matrix by Individual IDs
  common_ids <- intersect(rownames(G), names(phenotype))

  if (length(common_ids) == 0) {
    stop("No common individual IDs found between genotype matrix row names and phenotype vector names. Cannot align data.")
  }

  if (length(common_ids) < nrow(G) || length(common_ids) < length(phenotype)) {
    warning(paste0("Only ", length(common_ids), " common individuals found. Some individuals in genotype or phenotype data do not have a match in the other dataset and will be excluded. Final number of individuals: ", length(common_ids)))
  }

  # Subset and order G and phenotype to include only common, aligned IDs
  G_aligned <- G[common_ids, common_ids]
  phenotype_aligned <- phenotype[common_ids]

  # Additional checks on aligned data
  if (length(phenotype_aligned) == 0) {
      stop("After alignment, no individuals remain. Please check your IDs and data.")
  }
  if (var(phenotype_aligned, na.rm=TRUE) == 0) {
      stop("Phenotype data has zero variance after alignment. Heritability cannot be calculated.")
  }
  if (any(is.nan(phenotype_aligned)) || any(is.infinite(phenotype_aligned))) {
      stop("Phenotype data contains NaN or Inf values after alignment. Please check your raw data for errors or inconsistencies.")
  }
  if (any(is.na(phenotype_aligned))) {
      message("Phenotype data contains NA values. BGLR will handle these during estimation by excluding them from the likelihood.")
  }


  # 4. Define the BGLR Model
  ETA <- list(
    list(K = G_aligned, model = "BRR")
  )

  # 5. Run BGLR
  message("Running BGLR model... This can take some time depending on nIter.")
  fit <- BGLR(
    y = phenotype_aligned,
    ETA = ETA,
    nIter = nIter,
    burnIn = burnIn,
    thin = thin,
    verbose = verbose,
    saveAt = "BGLR_output_" # Optional: prefix for saving outputs if needed
  )

  # 6. Extract Variance Components and Perform Diagnostics
  var_additive_samples <- fit$ETA[[1]]$varB
  var_residual_samples <- fit$varE

  # --- DIAGNOSTIC CHECKS ON BGLR OUTPUT SAMPLES ---
  # Check if samples are empty or contain problematic values
  if (length(var_additive_samples) == 0 || any(is.na(var_additive_samples)) || any(is.infinite(var_additive_samples))) {
    warning("Additive variance (varB) samples are empty or contain NA/Inf values after BGLR estimation. This suggests issues with model convergence or data. Returning NA for heritability.")
    var_additive <- NA
  } else {
    var_additive <- mean(var_additive_samples)
  }

  if (length(var_residual_samples) == 0 || any(is.na(var_residual_samples)) || any(is.infinite(var_residual_samples))) {
    warning("Residual variance (varE) samples are empty or contain NA/Inf values after BGLR estimation. This suggests issues with model convergence or data. Returning NA for heritability.")
    var_residual <- NA
  } else {
    var_residual <- mean(var_residual_samples)
  }

  # 7. Calculate Narrow-Sense Heritability
  # Ensure variance components are valid numbers before calculation
  if (is.na(var_additive) || is.na(var_residual) || (var_additive + var_residual) == 0) {
      h2_narrow <- NA # Cannot calculate heritability if variances are problematic
      warning("Cannot calculate h2: Variance components are NA/problematic or sum to zero after BGLR estimation.")
  } else {
      h2_narrow <- var_additive / (var_additive + var_residual)
  }

  # 8. Return Results
  results <- list(
    h2_narrow = h2_narrow,
    var_additive = var_additive,
    var_residual = var_residual,
    BGLR_model = fit # Return the full BGLR model object for detailed inspection
  )

  message("Calculation complete.")
  return(results)
}

# Helper function to check for positive semi-definiteness
is.positive.semidefinite <- function(m) {
  if (!is.matrix(m) || nrow(m) != ncol(m)) return(FALSE)
  eigen_values <- tryCatch(eigen(m, symmetric = TRUE, only.values = TRUE)$values,
                           error = function(e) {
                             warning("Eigenvalue calculation error for G matrix during PSD check: ", e$message)
                             return(NULL)
                           })
  if (is.null(eigen_values)) {
    return(FALSE)
  }
  return(all(eigen_values >= -1e-9))
}
#----------------------------------------------------------------
#-- Calculata heritability for a single trait using BayesA method
#----------------------------------------------------------------
new_calculateHeritabilityTrait <- function(traitName, geno, pheno, outputDir) {
	message(">>> Calculating Heritability for ", traitName)

	# Prepare data
	X = scale(geno) / sqrt(ncol(geno))	# Standardize genotypes
	Y = pheno[, traitName]				# Extract phenotype vector

	# Set model and output prefix
	MODEL = "RKHS"
	PREFIX = paste0(outputDir, "/tmp/", MODEL, "-", traitName, "-")

	# Create kernel matrix for genetic relationship
	K = tcrossprod(X) / ncol(X)  # Centered GRM (Genomic Relationship Matrix)

	# Fit the model using BGLR
	eta = list(list(K=K, model=MODEL))
	fm = BGLR(y=Y, ETA=eta, nIter=5000, burnIn=1000, verbose=FALSE, saveAt=PREFIX)

	# Extract variance components directly from the model object
	varU = fm$ETA[[1]]$varU  # Genetic variance
	varE = fm$varE			 # Residual variance

	# Calculate narrow-sense heritability
	print (paste0 ("VarU:", varU))
	print (paste0 ("VarE:", varE))
	h2 = varU / (varU + varE)

	# Return heritability estimate
	return(h2)
}

calculateHeritabilityTrait <- function (traitName, geno, pheno, outputDir) {
	message (">>> Heritability for ", traitName)
	# Get Trait and Geno for BGLR
	X = scale(geno)/sqrt(ncol(geno))
	Y = pheno [, traitName]
	 
	# Estimation using variance components
	MODEL  = "BRR"
	PREFIX = paste0 (outputDir, "/tmp/", MODEL, "-", traitName,"-")
	varUFilename = paste0 (PREFIX, "ETA_1_varB.dat")
	#varUFilename = paste0 (PREFIX, "ETA_1_ScaleBayesA.dat")
	varEFilename = paste0 (PREFIX, "varE.dat")

	h2_0 = .5
	eta  = list(list(X=X,model=MODEL,saveEffects=T))
	fm	 = BGLR(y=Y,ETA=eta, nIter=5000, burnIn=1000, verbose=F, saveAt=PREFIX)

	varU = scan(varUFilename)
	varE = scan(varEFilename)

	# Calculate narrow-sense heritability
	print (paste0 ("VarU:", varU))
	print (paste0 ("VarE:", varE))

	h2	 = varU/(varU+varE)

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
		Prefix		  = strsplit (traitNames[i], "[.]")[[1]][1]
		HCL_component = strsplit (traitNames[i], "[.]")[[1]][2]
		Models		  = "BA"
		Traits		  = traitNames [i]
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

