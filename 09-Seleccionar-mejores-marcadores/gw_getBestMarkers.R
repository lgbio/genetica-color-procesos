#!/usr/bin/Rscript
library (dplyr)
library (ggmanh)
library (ggrepel) # For SNP labels
source ("lglib14.R")
source ("gw_scoreMarkers.R")
warnings()
	
#-------------------------------------------------------------
# Best markers for both base and HCL traits
# For each trait N best markers are selected
# Create manhattan plots for top best markers and 
#-------------------------------------------------------------

NBEST = 4634          # 4634 SNPs
INPUT = "inputs/SNPs-GSCORES-ALL-3TOOLS.csv" # GWAS results for 4634 SNPs, 9 BASE TRAITS, 3 HCL TRAITS, 3 TOOLS 
OUT_ALL_BASE  = "outputs/SNPs-GSCORES-3TOOLS-BEST-ALL-SNPs-BASE.csv" # Best 4634 SNPs for each base trait (42156 SNPs)
OUT_ALL_HCL = "outputs/SNPs-GSCORES-3TOOLS-BEST-ALL-SNPs-HCL.csv"  # Best 3 HCL SNPs for each HCL trait (42156*3 SNPs)

OUT_TOP_BASE  = "outputs/SNPs-GSCORES-3TOOLS-BEST-TOP-SNPs-BASE.csv" # Best 4634 SNPs for each base trait (42156 SNPs)
OUT_TOP_HCL = "outputs/SNPs-GSCORES-3TOOLS-BEST-TOP-SNPs-HCL.csv"  # Best 3 HCL SNPs for each HCL trait (42156*3 SNPs)

main <- function () {
    doAnalysis (INPUT)
}

doAnalysis <- function (inputData) {
    createDir ("outputs")
#reload (pkgload::inst("gwascolors"))
	markersFile = inputData

	# Add column with base trait name
	markersFile = addBaseTraits (markersFile)                     #

    # Best N for base traits
	outFile = gw_bestMarkersTraits (markersFile, NBEST, "BASE")
	file.symlink (basename (outFile), OUT_ALL_BASE)

	# Best N for HCL traits
	outFile = gw_bestMarkersTraits (markersFile, NBEST, "TRAIT")
	file.symlink (basename (outFile), OUT_ALL_HCL)

    # Best one for plotting
	outFile = gw_bestMarkersTraits (markersFile, 1, "BASE")
	file.symlink (basename (outFile), OUT_TOP_BASE)
	plotBestSNPsByTrait (outFile, "BASE")

    # Best one for Plotting
	outFile = gw_bestMarkersTraits (markersFile, 1, "TRAIT")
	file.symlink (basename (outFile), OUT_TOP_HCL)
	plotBestSNPsByTrait (outFile, "TRAIT")

	
#	# Best 10 SNPs for 
	outFile = gw_bestMarkersTraits (markersFile, 10, "BASE")
	plotBestSNPsByTrait (outFile, "BASE")
	outFile = gw_bestMarkersTraits (markersFile, 10, "TRAIT")
	plotBestSNPsByTrait (outFile, "TRAIT")
    #-----------------------------------------------------------------------------

	# Manhattan plot for top 1000 and labels for 20
	topFile        = gw_getTopNBestMarkers (markersFile, 1000)
	gw_plotManhattanBestMarkers (topFile, N=20)

	#gw_plotManhattanBestMarkers  (markersFile, N=20) 
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#'
#' Select best GWAS markers
#'
#' Select best N GWAS markers for each trait from GWAS results file. Selection is done based on a own score called "AgroScore" that takes into account three terms: Genomic Control (80%), replicability (10%) and significance (10%). 
#' Best markers are selected by trait: first, each of the four tools selects its N=50 best markers, low p-value, resulting in M=200 markers for each phenotype. These markers are scored with our "AgroScore" and the best N=50 markers are selected.
#'
#' @param  markersFile Filename of GWAS results.
#' @param  nBest       Numbers of best SNPs for each trait.
#' @param  TRAITTYPE   Type of trait to get best markers.
#' @param  outFile     Filename of resulting output file.
#' @return Table with best N markers for each trait. 
#' @export
#-------------------------------------------------------------
gw_bestMarkersTraits <- function (markersFile, nBest, TRAITTYPE) {
    nSuffix = sprintf ("-BEST-N%.2d-%s.csv", nBest, TRAITTYPE)
	outFile = gsub (".csv", nSuffix, markersFile)

	message ("Getting best N=", nBest, " markers for ", TRAITTYPE, " traits...")
	markers    = read.csv (markersFile); #view (markers,m=0)

	traitNames = unique (markers[,c(TRAITTYPE)])
	bestMarkersTable = NULL
#	for (trait in traitNames) {
#		scoresTrait = filter (markers, !!sym(TRAITTYPE)==trait);#view (scoresTrait)
#		uniqueMarkers  = distinct (scoresTrait, SNP, .keep_all=T);view (uniqueMarkers)
#		bestMarkers    = filter (uniqueMarkers, SIGNIFICANCE==T & GSCORE > GVALUE)
#		bestMarkersTable = rbind (bestMarkersTable, bestMarkers)
#	}
	top_snps <- markers %>%
	  # Group by Class (and Subclass if needed)
	  group_by (!!sym(TRAITTYPE)) %>%
	  # Sort by Score (descending) to prioritize high scores
	  arrange(desc(GSCORE), .by_group = TRUE) %>%
	  # Keep only distinct SNPs (no duplicates per group)
	  distinct (SNP, .keep_all = TRUE) %>%
	  # Slice the top N rows per group
	  slice_head (n = nBest) %>%
	  # Remove grouping (optional)
	  ungroup()

#    bestMarkersTable = markers %>% 
#        group_by (!!sym(TRAITTYPE)) %>% 
#        slice_max (GSCORE, n=nBest, with_ties=F) %>% 
#        ungroup ()

    outFile = gsub ("inputs", "outputs", outFile)
	write.csv (top_snps, outFile, row.names=F)
	return (outFile)
}
#-------------------------------------------------------------
#' Select best GWAS markers from all traits
#' @param  markersFile  Filename of GWAS results.
#' @param  outFile     Filename of resulting output file.
#' @return Output filename.
#' @export
#-------------------------------------------------------------
gw_getTopNBestMarkers <- function (markersFile, N=50, outFile="") {
	if (outFile=="") 
			outFile = gsub (".csv", sprintf ("-Top%s.csv", N), markersFile)

	markers       = read.csv (markersFile); #view (markers,m=0)
	uniqueMarkers = distinct (markers, SNP, .keep_all=T)
	if (nrow (uniqueMarkers) > N)
		uniqueMarkers = uniqueMarkers[1:N,]

    outFile = gsub ("inputs", "outputs", outFile)
	write.csv (uniqueMarkers, outFile, row.names=F)
	return (outFile)
}

#-------------------------------------------------------------
#-- Create Manhattan plot and label the first N SNPs
#-------------------------------------------------------------
gw_plotManhattanBestMarkers <- function (markersFile, N=20) {
	outFile     = gsub (".csv", "-Manhattan.pdf", markersFile)
	outTitle    = paste0 ("Top ", N, " best markers")
	markers     = read.csv (markersFile)

	# Create vector of N names for manhattan_plot
	names       = rep ("", nrow (markers))
	names [1:N] = markers$SNP [1:N]
	markers     = cbind (NAMES=names, markers);#view (markers)
	manhattan_plot (x=markers, pval.colname="P", chr.colname="CHR", pos.colname="POS", 
					plot.title=outTitle, y.label="Score\n-log10(P)", label.colname="NAMES",
					label.font.size=4, point.size=2, signif = c(1,1))
    outFile = gsub ("inputs", "outputs", outFile)
	ggsave (outFile, width=10, height=7)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
plotTraitsMultipleManhattanPlots <- function (markersFile, N=20) {
	library (cowplot)
	outFile  	  = gsub (".csv", "-Plot.pdf", markersFile)
	outTitle    = paste0 ("Top ", N, " best markers")
	markers     = read.csv (markersFile)

	# Write to file Top N markers
	outFileN    = gsub (".csv", "-N20.csv", markersFile)
    outFileN = gsub ("inputs", "outputs", outFileN)
	write.csv (markers [1:N,], outFileN, row.names=F)

	traits = unique (markers$TRAIT)
	gList = list()	
	for (t in traits[1:7]) {
			markersTrait = markers %>% filter (TRAIT==t)
			# Create vector of N names for manhattan_plot
			names       = rep ("", nrow (markersTrait))
			names [1:N] = markersTrait$SNP [1:N]
			markersTrait     = cbind (NAMES=names, markersTrait);#view (markersTrait)

			g = manhattan_plot (x=markersTrait, pval.colname="P", chr.colname="CHR", 
								pos.colname="POS", y.label="Score\n-log10(P)", 
								signif = c(1,1))
			gList <- append (gList, list(g))
	}
	plot_grid (plotlist=gList, ncol=1, axis="l", align="v")
    outFile = gsub ("inputs", "outputs", outFile)
	ggsave (outFile)
}
#-------------------------------------------------------------
# Plot best SNPs by trait in a Traits vs Chromosomes plot 
#-------------------------------------------------------------
plotBestSNPsByTrait <- function (markersFile, TRAITTYPE) {
	message ("Plotting best markers for ", TRAITTYPE, " traits...")
	yLabel = ifelse (TRAITTYPE=="BASE", "BASE TRAITS", "HCL TRAITS")
	outFile = gsub (".csv", ".pdf", markersFile)

	markers = read.csv (markersFile) %>% rename ("BP"="POS") 

	# First of all, we need to compute the cumulative position of SNP.
	markersBP <- markers %>% 
	  # Compute chromosome size
	  group_by(CHR) %>% summarise(chr_len=max(BP)) %>% 
	  # Calculate cumulative position of each chromosome
	  mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>%
	  # Add this info to the initial dataset
	  left_join(markers, ., by=c("CHR"="CHR")) %>%
	  # Add a cumulative position of each SNP
	  arrange(CHR, BP) %>% mutate( BPcum=BP+tot)

	# Then we need to prepare the X axis. Indeed we do not want to display 
	# the cumulative position of SNP in bp, but just show the chromosome name instead.
	axisdf = markersBP %>% group_by(CHR) %>% summarize(center=(max(BPcum)+min(BPcum))/2)

	# Add labels only to markers that are significant
	markersBP = markersBP %>% mutate (names=ifelse (SIGNIFICANCE, TRUE, FALSE))

	# Ready to make the plot using ggplot2:
	ggplot(markersBP, aes(x=BPcum, y=!!sym(TRAITTYPE))) +
	    # Show all points
	    geom_point( aes(color=as.factor(CHR)), alpha=0.6, size=1.5) +
	    scale_color_manual(values = rep(c("orange", "midnightblue"), 22 )) +
	    # custom X axis:
	    scale_x_continuous (label = axisdf$CHR, breaks= axisdf$center ) +
	    # SNP labels
	    geom_label_repel (data=subset (markersBP, names==TRUE), 
	    				  aes (fill=factor(SNP), alpha=0.0, label=SNP), size=2) +
	    # Custom the theme:
	    theme_bw() + theme(legend.position="none",
	      panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
	      panel.border = element_rect(colour = "black", fill=NA)) +
        labs (title="Best SNP by trait", y=yLabel, x="CHROMOSOME")
    outFile = gsub ("inputs", "outputs", outFile)
	ggsave (outFile, width=7, height=4)
}

#-------------------------------------------------------------
# Add base trait from HCL trait
#-------------------------------------------------------------
addBaseTraits <- function (markersFile) {
	message ("Adding base trait...", markersFile)
	outFile = gsub (".csv", "-TRAITS.csv", markersFile)
	markers    = read.csv (markersFile); #view (markers,m=0)
	baseTraits = sapply (markers$TRAIT, function (x) strsplit (x, "[.]")[[1]][1])
	markers    = cbind (BASE=baseTraits, markers)

    outFile = gsub ("inputs", "outputs", outFile)
	write.csv (markers, outFile, row.names=F)
	return (outFile)
}



#-------------------------------------------------------------
#-------------------------------------------------------------
main()
