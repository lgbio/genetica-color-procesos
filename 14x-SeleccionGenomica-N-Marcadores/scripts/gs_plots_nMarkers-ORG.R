#!/usr/bin/env Rscript
source ("lglib14.R")   # LuisG library with utility functions

OUTFILETABLE = "outputs/out-GS-crossval-varying-nmarkers.csv"
OUTBASE      = "outputs/out-GS-crossval-varying-nmarkers-basetraits.pdf"
OUTHCLs      = "outputs/out-GS-crossval-varying-nmarkers-hcltraits.pdf"

#-------------------------------------------------------------
# Function for ploting predictive ability varying number of markers
#-------------------------------------------------------------
main <- function () {
	    #library (ppGS)
	library (ggplot2)
	library (dplyr)

    createDir ("outputs")
	args = c ("runs")

	runsDir = args [1]

	summaryTableFile = gs_plots_create_summary_table (runsDir, OUTFILETABLE)
	system ("change-traitnames.py")
	gs_plots_nmarkers_summary (summaryTableFile, OUTBASE, OUTHCLs)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#' Create summaries when the number of markers varies. 
#'
#' Create one summary table and two plots: one per base trait and one per hcl trait.
#' @param inputDir  Dir name with results from different GS runs.
#' @return None. Writes table (.csv) and plots (.pdfs) to files.
#' @import dplyr ggplot2
#' @export
gs_plots_nmarkers_summary <- function (summaryFile, outBaseFile, outHCLsFile) {
	summaryTableAll      = read.csv (summaryFile)
	summaryTableFiltered = filter_traits (summaryTableAll)
	view (summaryTableFiltered)

	# Define X breaks
	runsVector = unique (summaryTableAll$nMarkers)
    runsBreaks = sort (as.integer (runsVector))
    runsLabels = paste0 (runsBreaks, "")
    runsLimits = c (0, max (runsBreaks))

	#-- First: Plot for base traits --------------------------------------------------------
	ggplot (summaryTableFiltered, aes(x=nMarkers, y=Means, group=HCLs)) + #ylim (0,1) +
			geom_line (aes (color=HCLs)) +
			geom_point (aes (color=HCLs)) +
			labs (title="Genomic selection with varying marker densities per trait") +
			ylab ("Predictive ability") + xlab ("Percentage of markers (%)") +
			scale_x_continuous(breaks=runsBreaks, labels=runsLabels, limits=runsLimits) +
			facet_wrap (~Traits, ncol=3) 

	ggsave (outBaseFile, width=11)

	#-- Second: Plot for base traits --------------------------------------------------------
	ggplot (summaryTableFiltered, aes(x=nMarkers, y=Means, group=Traits)) + 
		#ylim (0,1) +
		geom_line (aes (color=Traits)) +
		geom_point (aes (color=Traits)) +
			labs (title="Genomic selection with varying marker densities per HCL trait") +
		ylab ("Predictive ability") + xlab ("Percentage of Markers (%)") +
		scale_x_continuous(breaks=runsBreaks, labels=runsLabels, limits=runsLimits) +
		facet_wrap (~HCLs, ncol=3) +
		theme_bw() +	# This changes to a white background with light gray gridlines
		theme(
			panel.background = element_rect(fill = "white"),  # Ensures pure white background
			plot.background = element_rect(fill = "white"),   # White background for entire plot
			panel.grid.major = element_line(color = "grey90"), # Lighter grid lines
			panel.grid.minor = element_line(color = "grey95")  # Even lighter minor grid lines
		)		   

		ggsave (outHCLsFile, width=10, height=7, bg='white')
}

#-------------------------------------------------------------
# Remove specific traits for publication
#-------------------------------------------------------------
filter_traits <- function (summaryTable) {
	# Define the %notin% operator
	`%notin%` <- Negate(`%in%`)
	# Traits to remove
	traits_to_remove <- c("BerryC", "PCTuberflesh")

	# Using %notin%
	df_filtered <- summaryTable [summaryTable$Traits %notin% traits_to_remove, ]
	return (df_filtered)
}


#-------------------------------------------------------------
##' Create summaries when the number of markers varies. 
#'
#' Create one summary table 
#' @param inputDir  Dir name with results from different GS runs.
#' @return None. Writes table (.csv) 
#' @import dplyr 
#-------------------------------------------------------------
gs_plots_create_summary_table <- function (inputDir, outFile) {
	#library (dplyr)
	#library (ggplot2)

	# Create summary table
	dirList = list.dirs (inputDir, recursive=F)
	summaryTableList = list()
    runsVector = c()
	for (runDir in dirList) {
        outDir = sprintf ('%s/%s', runDir, 'outputs')
		print (outDir)
		#n = as.numeric (gsub (paste0(inputDir,"/out"), "", outDir))
		n = as.numeric (unlist (strsplit (outDir, "/"))[2])
        runsVector = c (runsVector, n)
		hclData = read.csv (sprintf ("%s/%s", outDir, "out-CrossValidation-kfolds-GEBVs.csv"))
		summaryTable = data.frame (nMarkers=n, createSummaryTableKFolds (hclData))
		summaryTableList = append (summaryTableList, list (summaryTable))
	}
	summaryTableAll = do.call ("rbind", summaryTableList)
	write.csv (summaryTableAll, outFile, row.names=F)
	#write.csv (summaryTableAll, "outputs/out-GS-crossval-varying-nmarkers.csv", row.names=F)
	return (outFile)
}

#-------------------------------------------------------------
#' Create summaries when the number of markers varies. 
#'
#' Create one summary table and two plots: one per base trait and one per hcl trait.
#' @param inputDir  Dir name with results from different GS runs.
#' @return None. Writes table (.csv) and plots (.pdfs) to files.
#' @import dplyr ggplot2
#' @export
old_gs_plots_varying_nmarkers <- function (inputDir) {
	#library (dplyr)
	#library (ggplot2)

	# Create summary table
	dirList = list.dirs (inputDir, recursive=F)
	summaryTableList = list()
    runsVector = c()
	for (runDir in dirList) {
        outDir = sprintf ('%s/%s', runDir, 'outputs')
		print (outDir)
		#n = as.numeric (gsub (paste0(inputDir,"/out"), "", outDir))
		n = as.numeric (unlist (strsplit (outDir, "/"))[2])
        runsVector = c (runsVector, n)
		hclData = read.csv (sprintf ("%s/%s", outDir, "out-CrossValidation-kfolds-GEBVs.csv"))
		summaryTable = data.frame (nMarkers=n, createSummaryTableKFolds (hclData))
		summaryTableList = append (summaryTableList, list (summaryTable))
	}
	summaryTableAll = do.call ("rbind", summaryTableList)
	write.csv (summaryTableAll, "outputs/out-GS-crossval-varying-nmarkers.csv", row.names=F)

	# Create plot per base trait
    runsBreaks = sort (as.integer (runsVector))
    print (runsBreaks)
    runsLabels = as.character (runsBreaks)
    print (runsLabels)
    runsLimits = c (0, max (runsBreaks))
    print (runsLimits)

	ggplot (summaryTableAll, aes(x=nMarkers, y=Means, group=HCLs)) + #ylim (0,1) +
			geom_line (aes (color=HCLs)) +
			geom_point (aes (color=HCLs)) +
			labs (title="Genomic selection varying the number of markers per base trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			#scale_x_discrete (name = "Number of Markers", limits = c ("5", "15", "25", "35", "50", "100")) +
			#scale_x_continuous(breaks = c(5,15,25,35,50,100), labels = c("5","15","25","35","50","100"), limits = c(0,100)) +
			scale_x_continuous(breaks=runsBreaks, labels=runsLabels, limits=runsLimits) +
			facet_wrap (~Traits, ncol=3) 

	ggsave ("outputs/out-GS-crossval-varying-nmarkers-basetraits.pdf", width=11)

	# Create plot per hcl trait
	ggplot (summaryTableAll, aes(x=nMarkers, y=Means, group=Traits)) + 
			#ylim (0,1) +
			geom_line (aes (color=Traits)) +
			geom_point (aes (color=Traits)) +
			labs (title="Genomic selection varying the number of markers per HCL trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			#scale_x_continuous(breaks = c(5,15,25,35,50,100), labels = c("5","15","25","35","50","100"), limits = c(0,100)) +
			scale_x_continuous(breaks=runsBreaks, labels=runsLabels, limits=runsLimits) +
			facet_wrap (~HCLs, ncol=3) 


	ggsave ("outputs/out-GS-crossval-varying-nmarkers-hcltraits.pdf", width=10, height=7)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#inputDir = "outs"
#library (dplyr)
#library (ggplot2)
#gs_plots_varying_nmarkers (inputDir)


#-------------------------------------------------------------
#--- Create output tables and boxplots for best model for each trait
#--- Return a summary table with hcl trait, baseTrait, hclComp, best model, and means
#-------------------------------------------------------------
createSummaryTableKFolds <- function (HCLTable) {
	message (">>> Creating summary table of kfolds: hclTrait, baseTrait, hclComp, model, means ...")

	# Summary plot boxplot best models for traits
	meansTable = HCLTable %>% dplyr::group_by (Traits, Models) %>% 
		summarize (Means=mean(get("Predictive_ability")), .groups="keep") 

	# Get best model for each trait
	bestTable = meansTable %>% group_by (Traits) %>% slice_max (Means) %>% slice_head

	# Create table with correlations from best models
	getBestCorr <- function (trait, model) 
		return (dplyr::filter (HCLTable, Traits==trait, Models==model))
	bestCorrs = do.call (rbind, mapply (getBestCorr, bestTable$Traits, bestTable$Models, SIMPLIFY=F))

	# Add base trait and HLC componente to summary table
	baseTraitList = c ()
	hclCompList   = c ()
	for (hclTrait in bestTable$Traits){
		baseTraitList = c (baseTraitList, strsplit (hclTrait, "[.]")[[1]][1])
		hclCompList   = c (hclCompList, strsplit (hclTrait, "[.]")[[1]][2])
	}
	summaryTable = data.frame (bestTable[,1], baseTraitList, hclCompList, bestTable[,2], bestTable [,3])
	names (summaryTable) = c("HCLTraits", "Traits", "HCLs", "Models", "Means")
	write.csv (summaryTable, "outputs/out-CrossValidation-GEBVs-Summary.csv", row.names=F)
	return (summaryTable)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
main ()
