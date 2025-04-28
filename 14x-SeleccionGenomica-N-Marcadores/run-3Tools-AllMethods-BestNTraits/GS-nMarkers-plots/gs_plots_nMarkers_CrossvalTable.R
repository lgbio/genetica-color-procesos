#!/usr/bin/env Rscript

#-------------------------------------------------------------
# Function for ploting predictive ability varying number of markers
# from GS crossvalidation table file
#-------------------------------------------------------------
INPUT	= "out-GS-crossval-varying-nmarkers.csv"
OUTBASE = "out-GS-crossval-varying-nmarkers-basetraits.pdf"
OUTHCLs = "out-GS-crossval-varying-nmarkers-hcltraits.pdf"

#-------------------------------------------------------------
#-------------------------------------------------------------
main <- function () {
	source ("lglib14.R")   # LuisG library with utility functions
	library (ggplot2)
	library (dplyr)

	gs_plots_nmarkers_summary (INPUT, OUTBASE, OUTHCLs)
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
#-------------------------------------------------------------
gs_plots_nmarkers_summary <- function (summaryFile, outBaseFile, outHCLsFile) {
	#library (dplyr)
	#library (ggplot2)

	# Read summary table
	summaryTableAll		 = read.csv (summaryFile, check.names = F)
	summaryTableFiltered = filter_traits (summaryTableAll)
	view (summaryTableFiltered)

	# Define plot breaks, plot labels, and plot limits
	runsVector = c(5, 10, 25, 50, 75, 100)	 # Runs for N Markers
	runsBreaks = sort (as.integer (runsVector))
	runsLabels = as.character (runsBreaks)
	runsLimits = c (0, max (runsBreaks))

	#-- First: Plot for base traits --------------------------------------------------------
	ggplot (summaryTableFiltered, aes(x=nMarkers, y=Means, group=HCLs)) + #ylim (0,1) +
			geom_line (aes (color=HCLs)) +
			geom_point (aes (color=HCLs)) +
			labs (title="Genomic selection varying the number of markers per base trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			scale_x_continuous(breaks=runsBreaks, labels=runsLabels, limits=runsLimits) +
			facet_wrap (~Traits, ncol=3) 

	ggsave (outBaseFile, width=11)

	#-- Second: Plot for base traits --------------------------------------------------------
	ggplot (summaryTableFiltered, aes(x=nMarkers, y=Means, group=Traits)) + 
		#ylim (0,1) +
		geom_line (aes (color=Traits)) +
		geom_point (aes (color=Traits)) +
		labs (title="Genomic selection varying the number of markers per HCL trait") +
		ylab ("Predictive ability") + xlab ("Number of markers") +
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
#-------------------------------------------------------------
main ()
