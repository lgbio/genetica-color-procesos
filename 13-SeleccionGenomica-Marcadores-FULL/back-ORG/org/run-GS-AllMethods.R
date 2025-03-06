#!/usr/bin/Rscript
detach("package:ppGS", unload=TRUE)
library (ppGS)
# Get sourde geno and pheno datasets (pheno: componentsLCH)
gs_datasets ("sources")

# Create configuration file multiple phenotypes
gs_configfile (type="multi")
set.seed (123)
gs_split ("sources/sources-genotypes-CCC-Andigena-ClusterCall2020.csv", 
		  "sources/sources-phenotypes-CCC-Andigena-ColorTraitsHCL.csv",
		  outputDir="data")
#, nGenos=1000, nPheno=3)
gs_multi ("config-multi.yml", "outputs")
