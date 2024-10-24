#!/usr/bin/Rscript
GOAL="
Get common SNPs by both base trait and derived trait.
Input markers are scored by AgroScore"

library (dplyr)
library (stringi)
source ("lglib14.R")
options (width=300)

`%notin%` <- Negate(`%in%`)
#-----------------------------------------------------------
#-----------------------------------------------------------
main <- function () {
    createDir ("outputs")
#	plotBestSNPsByTrait (outFile, labels=F)

    INPUTFILE        = "inputs/test.csv"
	#scores = read.csv ("inputs/SNPs-BESTN-TOOLS-ADDITIVE.csv", check.names=F)
	scores = read.csv ("inputs/SNPs-GSCORES-ALL-ADDITIVE-Top1000.csv", check.names=F)
	#scores = read.csv ("inputs/SNPs-GSCORES-ALL-ADDITIVE-TRAITS-Top1000.csv", check.names=F)
	#scores           = read.csv (INPUTFILE, check.names=F)
	#scores           = removeDuplicatedSNPsAllPhenos (scores);# view (scoresTable)
	commonBaseTraits = getCommonMarkersByBaseTrait (scores); #view (commonBaseTraits)
	commonHCLTraits  = getCommonMarkersByHCLTrait (scores); #view (commonHCLTraits)

}
#-----------------------------------------------------------
#-----------------------------------------------------------
getCommonMarkersByBaseTrait <- function (markers) {
	outFileWide = "outputs/common-SNPs-table-BaseTraits.csv"
#   outFileLong = "outputs/common-SNPs-table-BaseTraits-Long.csv"

	# Local function to join trait names
	xjoin  <- function (x) {stri_join (unique(x),collapse=",")}
	xcount <- function (x) {length (unique(x))}
	
	# Add base trait (BTRAIT)
	markers = markers %>% mutate (BTRAIT=sapply (strsplit (markers$TRAIT, "[.]"), function (x) x[1])) 

	# Select main columns and group by SNP
	markersSel   = markers %>% select (GSCORE, SNP, CHR, POS, BTRAIT, TRAIT)
	markersGrp   = markersSel %>% group_by (SNP) %>% distinct (SNP, BTRAIT, .keep_all=T)
	markersCount = markersGrp %>% mutate (nTRAITS=xcount (BTRAIT)) #%>% filter (nTRAITS > 1)

	# Results in long format 
#	write.csv (markersCount, outFileLong , row.names=F)

	# Results in wide format adding column with common traits by SNP
	markersWide = markersCount %>% mutate (COMMONTRAITS=xjoin(BTRAIT)) 
	markersWide = markersWide %>% distinct (SNP, .keep_all=T)
	markersWide$BTRAIT = NULL
	write.csv (markersWide, outFileWide, row.names=F)
	
	return (outFileWide)
}

#-----------------------------------------------------------
#-----------------------------------------------------------

plotBestSNPsByTrait <- function (markersFile, outFile="", labels=F) {
	library (ggrepel) # For SNP labels
	markers = read.csv (markersFile)

	if (outFile=="")
		outFile = gsub (".csv", ".pdf", markersFile)

	gwasResults = markers %>% rename ("BP"="POS") 

	# First of all, we need to compute the cumulative position of SNP.
	MAXBP = max (gwasResults$BP)
	don <- gwasResults %>% 
	  # Compute chromosome size
	  group_by(CHR) %>% summarise(chr_len=max(BP+BP/5)) %>% 
	  # Calculate cumulative position of each chromosome
	  mutate(tot=cumsum(chr_len)-chr_len) %>%
	  select(-chr_len) %>%
	  # Add this info to the initial dataset
	  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
	  # Add a cumulative position of each SNP
	  arrange(CHR, BP) %>% mutate( BPcum=BP+tot)
	
	# Then we need to prepare the X axis. Indeed we do not want to display 
	# the cumulative position of SNP in bp, but just show the chromosome name instead.
	axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

	# Add labels
	don = don %>% mutate (names=ifelse (labels,"yes", "no"))

	# Ready to make the plot using ggplot2:
	ggplot(don, aes(x=BPcum, y=BTRAIT)) +
	    # Show all points
	    #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.5) +
	    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.5) +
	    scale_color_manual(values = seq(1:56), 22 ) +
	    #scale_color_manual(values = rep(c("orange", "midnightblue"), 22 )) +
	    # custom X axis:
	    scale_x_continuous (label = axisdf$CHR, breaks= axisdf$center ) +
	    # SNP labels
	    geom_label_repel (data=subset (don, names=="yes"), aes (label=SNP), size=2) +
	    # Custom the theme:
	    theme_bw() + theme(legend.position="none",
	      panel.grid.major.x = element_blank(), 
	      panel.grid.minor.x = element_line (color="navyblue", linetype=2), 
	      panel.border = element_rect(colour = "black", fill=NA)) +
          labs (title="Common SNPs by base trait", y="BASE TRAITS", x="CHROMOSOMES")
	ggsave (outFile, width=7, height=7)
}
#-----------------------------------------------------------
#-----------------------------------------------------------
# Remove duplicated SNPs detected by the four GWAS toools
# Score each SNP and remove duplicates with low score 
#-----------------------------------------------------------
removeDuplicatedSNPsAllPhenos <- function (scores) {
	#--- Function for single phenos ---
	removeDuplitactedSNPsPheno <- function (gwasMarkers, phenoName) {
		scoresTable = gwasMarkers [gwasMarkers[,1]==phenoName,]

		# Score SNPs
		# Remove duplicated SNPs by keeping the SNP with highest score
		dups = duplicated (scores$SNP)
		snps = unique (scores [dups, "SNP"])
		scoresDups       = scores [scores$SNP %in% snps,]; #view (scoresDups)
		scoresDupsBest   = scoresDups %>% 
			               arrange (-GSCORE) %>% distinct (SNP, .keep_all=T); 
		scoresNoDups     = scores %>% filter (SNP %notin% snps); #view (scoresNoDups)
		scoresBest       = rbind (scoresDupsBest, scoresNoDups); #view (scoresBest)
		scoresBestSorted = scoresBest %>% arrange (TOOL, -GSCORE); #view (scoresBestSorted)
		return (scoresBestSorted)
	}
   
	phenosHCL      = unique (scores$TRAIT)
	scoresFiltered = lapply (phenosHCL, function (x) removeDuplitactedSNPsPheno (scores, x)) 
	scoresTable    = do.call ("rbind", scoresFiltered)
	return (scoresTable)
}

#-----------------------------------------------------------
# Get common markers for HCL traits (e.g. "Tallo.LCH.L", "Baya.LCH.H", ...)
#-----------------------------------------------------------
getCommonMarkersByHCLTrait <- function (scoresTable) {
	snpList    = unique (scoresTable$SNP)

	traitsSNPList = list()
	traitsSNPList2 = list()
	for (snp in snpList) {
		traits = scoresTable %>% filter (SNP==snp) %>% select (TRAIT)
		n  = nrow (traits)
		df = data.frame (SNP=snp, TRAIT=traits[,1], N=n) 
		traitsSNPList = append (traitsSNPList, list(df))

		df2 = data.frame (SNP=snp, TRAITS=paste (traits[,1], collapse=", "), N=n)
		traitsSNPList2 = append (traitsSNPList2, list(df2))
	}
	traitsSNPTable  = do.call ("rbind", traitsSNPList); #view (traitsSNPTable)
	traitsSNPSorted = traitsSNPTable %>% arrange (-N) ; #view (traitsSNPSorted)

	traitsSNPTable2  = do.call ("rbind", traitsSNPList2);# view (traitsSNPTable)
	traitsSNPSorted2 = traitsSNPTable2 %>% arrange (-N) ;# view (traitsSNPSorted)
	write.csv (traitsSNPSorted2, "outputs/common-SNPs-table-HCLTrait.csv", row.names=F)
	return (traitsSNPSorted)
}

#-------------------------------------------------------------------------
#-- Get common markers by base trait (e.g "Tallo", "Baya", "PriFlor", ...)
#-------------------------------------------------------------------------
old_getCommonMarkersByBaseTrait <- function (scoresTable) {
	allPhenosNames = sapply (strsplit (as.character (scoresTable$TRAIT), "[.]"), function (x) x[[1]][1])
	scoresPhenos   = data.frame (TRAITNAME=allPhenosNames, scoresTable);# view (scoresPhenos)

	snpList        = unique (scoresTable$SNP)

	traitsSNPList = list()
	traitsSNPList2 = list()
	for (snp in snpList) {
		traits = scoresPhenos %>% filter (SNP==snp) %>% select (TRAITNAME)
		traitsUnique = unique (traits [,1])
		n = length (traitsUnique)
		df  = data.frame (SNP=snp, TRAITS=traitsUnique, N=n) 
		traitsSNPList = append (traitsSNPList, list(df))

		df2 = data.frame (SNP=snp, TRAITS=paste (traitsUnique, collapse=", "), N=n)
		traitsSNPList2 = append (traitsSNPList2, list(df2))
	}
	traitsSNPTable  = do.call ("rbind", traitsSNPList);# view (traitsSNPTable)
	traitsSNPSorted = traitsSNPTable %>% arrange (-N) ;# view (traitsSNPSorted)

	traitsSNPTable2  = do.call ("rbind", traitsSNPList2);# view (traitsSNPTable)
	traitsSNPSorted2 = traitsSNPTable2 %>% arrange (-N) ;# view (traitsSNPSorted)
	write.csv (traitsSNPSorted2, "common-SNPs-table-BaseTrait.csv", row.names=F)
	return (traitsSNPSorted)
}

main()

