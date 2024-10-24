#!/usr/bin/env Rscript

# Create two tables: one with two codes for correlations. And
# the other with the final HCL traits and values for analizes
source ("lglib14.R")

USAGE="
Create table with the two color codifications for phenotypes from original colors table"

INPUTFILE    = "inputs/datos-AndigenaCCC-FenotiposColor-EscalaRHS-MERGED.csv"
INPUTCOLUMNS = "inputs/trait-color-ColumNames.csv"

OUTPUTFILE   = "outputs/out-colors-TwoCodes-Phenos.csv"
OUTPUTHCL    = "outputs/out-colors-RHSValues-BaseTraits.csv"
createDir ("outputs")


# Inputs/Outputs
args = commandArgs (trailingOnly=T)
dataAll = read.csv (INPUTFILE, na.strings=c(".","N"," ","","..", "0", "00"), check.names=F)  

# 1. Select columns from original data. Get columns with pairs of color codes: simple and RHS codes
dataHCLTraits = dataAll
names (dataHCLTraits) = dataAll [1,]; #view (dataHCLTraits)
columns = read.csv (INPUTCOLUMNS)

originalColumnNames = c("Registro", unlist(mapply(c, columns[,5], columns[,4], SIMPLIFY = FALSE)))
data  = dataHCLTraits [-1, originalColumnNames]; #view (data)

newColumnNames = c("Registro", unlist(mapply(c, columns[,3], columns[,2], SIMPLIFY = FALSE)))
names (data) = newColumnNames ; #view (data)

# 2. Format columns (create codes) and create two colors table for selected traits
ncols   = ncol (data)
i       = 2
colsHCL = c (newColumnNames[1])
while (i < ncols) {
	colsHCL = c (colsHCL, newColumnNames [i+1])
	data [,i]   = sapply (data[,i], function (x) strsplit (x, " ")[[1]][1], USE.NAMES=F)
	data [,i+1] = gsub (" ", "", data[,i+1], fixed=T)
	i = i+2
}

# Write table with two codes for correlations
write.csv (data, OUTPUTFILE, row.names=F)

# Write table with HCL values for the following analyzes
tableHCL = data [, colsHCL]
write.csv (tableHCL, OUTPUTHCL, row.names=F)

