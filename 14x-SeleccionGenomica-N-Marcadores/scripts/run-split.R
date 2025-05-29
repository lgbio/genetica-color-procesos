#!/usr/bin/env Rscript
#detach("package:ppGS", unload=TRUE)

Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa/14x-SeleccionGenomica-N-Marcadores")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
#library (ppGS)
devtools::load_all(file.path(Sys.getenv("RENV_PROJECT"), "lib/package-ppGS"))



# Get sourde geno and pheno datasets (pheno: componentsLCH)
#gs_datasets ("sources")

# Create configuration file multiple phenotypes
#gs_configfile (type="multi")
set.seed (123)
gs_split ("inputs/genotipo.csv",
		  "inputs/fenotipos.csv",
		  outputDir="datasets")
file.copy ("inputs/best-gwas-markers.csv", "datasets")
