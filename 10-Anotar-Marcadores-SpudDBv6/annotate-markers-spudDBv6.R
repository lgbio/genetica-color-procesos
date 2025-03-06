#!/usr/bin/env Rscript 

# Annotations of SNP markers resulting from GWAS or GS with SPUDDB Genome assembly (v6.1) annotations.
# "Improved genome assembly and annotation (v6.1) for the doubled monoploid potato DM 1-3 516 R44"
# at http://spuddb.uga.edu/dm_v6_1_download.shtml
#
# Input files:
#	1. File with markers (SNPs) to annotate: short SNP id instead large id  (e.g c1_7694 
#	   instead solcap_snp_c1_7694)
#	   Example file: "markers-scores.csv"
#	2. File with markers map: Chrommosome and Position
#	   Genome v6.1 file: "solcap_69k_SNPs_DM_v6_1_pos.txt"
#	3. File with gene annotations: Gen id with starting and ending position
#	   Genome v6.1 file: "DM_1-3_516_R44_potato.v6.1.hc_gene_models.gff3"
#	4. File with GO term assignments: Gen id with GO terms
#	   Genome v6.1 file: "DM_1-3_516_R44_potato.v6.1.working_models.iprscan_go_terms.txt"
#	5. File with InterPro assignments: Gen id with annotations
#	   Genome v6.1 file: "DM_1-3_516_R44_potato.v6.1.working_models.iprscan_output.tsv"
#
# Output file:
#	1. File with columns "SNP","REF","SNP_POS","SEQID","GOTERMS","SOURCE","ANNOT"

source ("lglib14.R")
library (stringi)
library (dplyr)
warnings()

INPUT_01 = "inputs/SNPs-GSCORES-ALL-3TOOLS-TRAITS-BASE.csv"
INPUT_02 = "inputs/SNPs-GSCORES-ALL-3TOOLS-TRAITS-HCL.csv"
createDir ("outputs")

main <- function () {
    doAnnotation (INPUT_01)
    doAnnotation (INPUT_02)
}

doAnnotation <- function (inputFile) {
    markersFile   = inputFile

    # SpudDb 6_1 inputs
    markersMap    = "annotations/markers-map.txt"
    genesFile     = "annotations/gene-annotations.gff3"
    gotermsFile   = "annotations/go-terms.txt"
    iprscanFile   = "annotations/iprscan-annotations.tsv"

    markersScores = read.csv (markersFile); #view (markersScores)

    N             = nrow (markersScores)
    outFile       = gsub (".csv", "-SPUDDB6_1.csv", markersFile)

    markersScores = markersScores %>% distinct (SNP, .keep_all=T)
    markersScores = markersScores [1:50,]
    markersMap    = read.csv (markersMap, sep="\t"); #view (markersMap)

    message (">>> Joining markersScores and markersMap...")
    markersScores = markersScores [, c ("GSCORE", "SNP", "TRAIT", "SIGNIFICANCE"), drop=F]
    markersMap = markersMap [!duplicated (markersMap$SNP_ID), c("SNP_ID", "REF", "SNP_POS")]
    markersMap$SNP_ID = gsub ("solcap_snp_", "", markersMap$SNP_ID);#view (markersMap)

    markers = merge (markersScores, markersMap, by.x="SNP", by.y="SNP_ID", all.x=T, sort=F); #view (markers)

    message (">>> Joining markers with genes...")
    library (sqldf)
    genes         = read.csv (genesFile, sep="\t", header=F); #view (genes)
    names (genes) = c ("SCAFFOLD", "MSU", "TYPE", "INI", "END", "NONE1", "STRAND", "NONE2", "INFO"); #view (genes)
    genesmRNA = genes [genes$TYPE=="mRNA",]; #view (genesmRNA)

    markersGenes = sqldf ("SELECT GSCORE, SNP, TRAIT, SIGNIFICANCE, REF, SNP_POS, INI, END, INFO
                           FROM markers
                           LEFT JOIN genesmRNA ON SNP_POS BETWEEN INI and END")
    SEQID = sapply  (strsplit (markersGenes$INFO, ";"), function (x) gsub ("ID=", "", x[[1]]))
    #view (SEQID)
    markersGenesIds = cbind (markersGenes, SEQID=SEQID); #view (markersGenesIds)

    # Add GO Terms
    message (">>> Joining markers with GO terms")
    goterms       = read.csv (gotermsFile, sep="\t", header=F); #view (goterms)
    names (goterms) = c ("SEQID", "GOTERMS")
    markersGenesIdsGOs = sqldf ("SELECT GSCORE, SNP, TRAIT, SIGNIFICANCE, REF, SNP_POS, INI, END, mk.SEQID, GOTERMS
                                 FROM markersGenesIds as mk
                                 LEFT JOIN goterms as gt
                                 ON mk.SEQID = gt.SEQID")
    #view (markersGenesIdsGOs)


    # Add Interpro annotations
    iprscan = read.csv (iprscanFile, header=F, sep="\t"); 
    names (iprscan) = c ("SEQID", "N1", "N2", "SOURCEANNOT", "SOURCE_ID", "ANNOT1", 
                         "N3","N4","N5","N6","N7","IPRID","ANNOT2", "GOTERMS")
    #view (iprscan)

    # Filter interpro annotations by Columns, no duplicated, and empty annotation
    iprscanSel = iprscan [, c("SEQID", "SOURCEANNOT", "ANNOT1", "ANNOT2")];#view (iprscanSel)
    iprscanNoDups = iprscanSel [!duplicated (iprscanSel),];#view (iprscanNoDups)
    iprscanFiltered = iprscanNoDups [iprscanNoDups$ANNOT1!="-" & iprscanNoDups$ANNOT2!="",]; #view (iprscanFiltered)
    ANNOT1 = iprscanFiltered$ANNOT1
    ANNOT2 = iprscanFiltered$ANNOT2
    ANNOT = ifelse (nchar (ANNOT1) > nchar (ANNOT2), ANNOT1, ANNOT2)  
    #view (ANNOT)
    iprscanFiltered$ANNOT1 = NULL
    iprscanFiltered$ANNOT2 = NULL
    iprscanAnnot = cbind (iprscanFiltered, ANNOT=ANNOT); #view (iprscanAnnot)
    iprscanAnnotNoDups = iprscanAnnot [!duplicated (iprscanAnnot),];#view (iprscanAnnotNoDups)

    message (">>> Joining markers with Interpro annotations...")
    markersIprscan = sqldf ("SELECT GSCORE, SNP, TRAIT, SIGNIFICANCE, REF, SNP_POS, mk.SEQID, mk.GOTERMS, SOURCEANNOT, ANNOT
                            FROM markersGenesIdsGOs as mk
                            LEFT JOIN iprscanAnnotNoDups as ip
                            ON mk.SEQID = ip.SEQID")
    #view (markersIprscan)

    write.csv (markersIprscan, "outputs/markers_annotations_SPUDDB_Genome6_1-DUPSNPS.csv", row.names=F)

    # Filter final annotations
    annot  = read.csv ("outputs/markers_annotations_SPUDDB_Genome6_1-DUPSNPS.csv")
    ## Remove NA annotations
    annotSNP = annot [!is.na (annot$SNP),]; #view (annotSNP, 10)

    ## Remove splicing sufix from SEQID
    seqidsf = strsplit (annotSNP$SEQID, fixed=F, "[.][0-9]$"); #view (seqidsf,10)
    annotSNP$SEQID = unlist (seqidsf);

    ## Select unique SEQID
    annotSeqid = annotSNP [!duplicated (annotSNP$SEQID),]; #view (annotSeqid,10)

    # Create table with one SNP marker and multiple values for each column
    getValues <- function (snp, ann, col) {
        values = ann [ann$SNP==snp, col]
        values = na.omit (values)
        values = unique (values)
        values = sprintf ('%s', stri_join (sprintf ("%s", values), collapse="|"))
    }

    snps     = unique (annotSeqid$SNP)
    snpsList = list()
    for (s in snps) {
        chrom    = getValues (s, annotSeqid, "REF")
        gscore   = getValues (s, annotSeqid, "GSCORE")
        signif   = getValues (s, annotSeqid, "SIGNIFICANCE")
        pos      = getValues (s, annotSeqid, "SNP_POS")
        seqids   = getValues (s, annotSeqid, "SEQID")
        traits   = getValues (s, annotSeqid, "TRAIT")
        goterms  = getValues (s, annotSeqid, "GOTERMS")
        annotsrc = getValues (s, annotSeqid, "SOURCEANNOT")
        annot    = getValues (s, annotSeqid, "ANNOT")
        snpsList = append (snpsList, list (list (gscore, s, chrom, pos, signif, seqids, traits,  goterms, annotsrc, annot)))
    }
    snpsdf = do.call (rbind.data.frame, snpsList)
    names (snpsdf) = c ("GSCORE", "SNP","CHROM", "POS","SIGNIFICANCE","SEQID","TRAIT","GOTERMS","SOURCEANNOT","ANNOT") 
    outFile = gsub ("inputs", "outputs", outFile)
    write.csv (snpsdf, outFile, row.names=F, quote=T)

    # Create summary table
    summCols = c("GSCORE", "SNP", "CHROM", "POS", "TRAIT", "ANNOT")
    summTable = snpsdf [1:N,summCols]
    selFirst <- function (x) {
        strsplit (x,"[|]")[[1]][1]
    }
    summTable [,"ANNOT"] = sapply (summTable [, c("ANNOT")], selFirst)
    summFile = gsub (".csv", "-SUMMARY.csv", outFile)
    summFile = gsub ("inputs", "outputs", summFile)
    write.csv (summTable, summFile, row.names=F)
}
#--------------------------------------------------
#--------------------------------------------------
main ()

