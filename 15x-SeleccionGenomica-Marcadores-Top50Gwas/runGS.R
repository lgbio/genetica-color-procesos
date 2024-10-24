#!/usr/bin/Rscript
# Run Genomic Selection using a config file and an output dir

library (ppGS)
set.seed (123)

args = commandArgs (trailingOnly=T)
configFile = args [1]
outDir     = args [2]
gs_multi (configFile, outDir)
