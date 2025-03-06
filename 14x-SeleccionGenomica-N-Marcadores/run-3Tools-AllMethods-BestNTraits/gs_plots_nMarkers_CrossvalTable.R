#!/usr/bin/env Rscript

#-------------------------------------------------------------
# Function for ploting predictive ability varying number of markers
# from GS crossvalidation table file
#-------------------------------------------------------------
INPUT   = "out-GS-crossval-varying-nmarkers.csv"
OUTBASE = "out-GS-crossval-varying-nmarkers-basetraits.pdf"
OUTHCLs = "out-GS-crossval-varying-nmarkers-hcltraits.pdf"

#-------------------------------------------------------------
#-------------------------------------------------------------
main <- function () {
	library (ggplot2)
	library (dplyr)

	gs_plots_nmarkers_summary (INPUT, OUTBASE, OUTHCLs)
}
#-------------------------------------------------------------
#-------------------------------------------------------------
gs_plots_nmarkers_summary <- function (summaryFile, outBaseFile, outHCLsFile) {
	#library (dplyr)
	#library (ggplot2)

	# Read summary table
    summaryTableAll = read.csv (summaryFile, check.names = F)

    # Define plot breaks, plot labels, and plot limits
    runsVector = c(5, 10, 25, 50, 75, 100)   # Runs for N Markers
    runsBreaks = sort (as.integer (runsVector))
    runsLabels = as.character (runsBreaks)
    runsLimits = c (0, max (runsBreaks))

    #-- First: Plot for base traits --------------------------------------------------------
	ggplot (summaryTableAll, aes(x=nMarkers, y=Means, group=HCLs)) + #ylim (0,1) +
			geom_line (aes (color=HCLs)) +
			geom_point (aes (color=HCLs)) +
			labs (title="Genomic selection varying the number of markers per base trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			scale_x_continuous(breaks=runsBreaks, labels=runsLabels, limits=runsLimits) +
			facet_wrap (~Traits, ncol=3) 

	ggsave (outBaseFile, width=11)

    #-- Second: Plot for base traits --------------------------------------------------------
	ggplot (summaryTableAll, aes(x=nMarkers, y=Means, group=Traits)) + 
			#ylim (0,1) +
			geom_line (aes (color=Traits)) +
			geom_point (aes (color=Traits)) +
			labs (title="Genomic selection varying the number of markers per HCL trait") +
			ylab ("Predictive ability") + xlab ("Number of markers") +
			scale_x_continuous(breaks=runsBreaks, labels=runsLabels, limits=runsLimits) +
			facet_wrap (~HCLs, ncol=3) 

	ggsave (outHCLsFile, width=10, height=7)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
main ()
