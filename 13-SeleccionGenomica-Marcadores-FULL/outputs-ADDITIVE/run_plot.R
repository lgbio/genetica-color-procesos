#!/usr/bin/env Rscript

#------------------- For calling from a renv environment -------------------------------------------------
Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa/14x-SeleccionGenomica-N-Marcadores")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
#---------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------
# Main only for testing
#-------------------------------------------------------------
main <-  function() {
    library (dplyr)
    library (ggplot2)
	source ("lglib14.R")
	source ("gs_plotsHCL_20.R")

	englishNamesFile = "traitnames-BASE.csv"
	outputDir		 = "outputs"
	hclTableFilename = "out-CrossValidation-kfolds-GEBVs.csv"

	hclTableFilename_SPANISH = hclTablePreprocessing (hclTableFilename, englishNamesFile, outputDir, "SPANISH")
    view (hclTableFilename_SPANISH)
	gs_plotsHCL (hclTableFilename_SPANISH, outputDir, "Componentes_CHL", "Habilidad_Predictiva", 
				 "Habilidad Predictiva para los 21 Componentes de Color CHL")

	hclTableFilename_ENGLISH = hclTablePreprocessing (hclTableFilename, englishNamesFile, outputDir, "ENGLISH")
	gs_plotsHCL (hclTableFilename_ENGLISH, outputDir, "CHL_Components", "Predictive_Ability", 
				 "Predictive Ability for the 21 CHL Color Componentes")
    #plotHCLs (newHclTableFilename, outputDir)
}

#--------------------------------------------------------------------------------
# Change to english names, remove two traits
#--------------------------------------------------------------------------------
hclTablePreprocessing <- function (hclTableFilename, englishNamesFile, outputDir, outType) {
	mappings <- read.csv (englishNamesFile) # Load the name mappings
  
	df <- read.csv (hclTableFilename)       # Load the target table

  	# Remove specific traits
	`%notin%` <- Negate(`%in%`)
	view (df)
	if (outType == "SPANISH") {
		df <- df %>% filter(Prefix %notin% c('CBaya', 'CPulpa'))
		names(df)[names(df) %in% c("HCL_component", "Predictive_ability")] <- c("Componentes_CHL", "Habilidad_Predictiva") # Rename columns "a" and "b" to "x" and "y"
	}else if (outType == "ENGLISH") {
		df <- df %>% filter(Prefix %notin% c('BerryC', 'PCTuberflesh'))
		names(df)[names(df) %in% c("HCL_component", "Predictive_ability")] <- c("CHL_Components", "Predictive_Ability") # Rename columns "a" and "b" to "x" and "y"
	}

	view (df)

	# Apply gsub replacements to all values
    for (i in seq_len(nrow(mappings))) {
      df[] <- lapply(df, function(col) {
        if (outType == "ENGLISH") 
            gsub( mappings$SpanishName[i], mappings$EnglishName[i], col, ignore.case = TRUE)
        else
            gsub( mappings$EnglishName[i], mappings$SpanishName[i], col, ignore.case = TRUE)
      })
    }  
 
	# Save the new table
	output_filename <- file.path (outputDir, addLabel (basename (hclTableFilename), outType))
	write.csv (df, output_filename , quote=F, row.names=F)
	return (output_filename)
}


old_hclTablePreprocessing <- function (hclTableFilename, englishNamesFile) {
  # Load the name mappings
  mappings <- read.csv (englishNamesFile)

  # Load the target table
  df <- read.csv (hclTableFilename)

  # Apply gsub replacements to all values
  for (i in seq_len(nrow(mappings))) {
    df[] <- lapply(df, function(col) {
      gsub(
        mappings$SpanishName[i],
        mappings$EnglishName[i],
        col,
        ignore.case = TRUE
      )
    })
  }  

  # Remove specific traits
  `%notin%` <- Negate(`%in%`)
  view (df)
  df <- df %>% filter(Prefix %notin% c('BerryC', 'PCTuberflesh'))
  view (df)

  # Save the new table
  output_filename <- gsub("TABLE.csv", "PREPRO-TABLE.csv", hclTableFilename)
  write.csv (df, output_filename , quote=F, row.names=F)
  return (output_filename)
}

#-------------------------------------------------------------
# Plot heritabilities dot, bar, and boxes
#-------------------------------------------------------------
plotHCLs <- function (hclTableFilename, outputDir) {
	gs_plotsHCL (hclTableFilename, outputDir, "Predictive_ability", 
				 "Genomic Prefiction for HCL traits")
}

#--------------------------------------------------------------------
#--------------------------------------------------------------------
main ()
