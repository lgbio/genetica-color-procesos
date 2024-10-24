#!/usr/bin/Rscript
#detach("package:ppGS", unload=TRUE)
library (ppGS)
# Get sourde geno and pheno datasets (pheno: componentsLCH)
#gs_datasets ("sources")

# Create configuration file multiple phenotypes
#gs_configfile (type="multi")
set.seed (123)
gs_split ("inouts/genotipo-AndigenaCCC-ClusterCall2020-MATRIX.csv",
		  "inouts/out-fenotiposColor-ComponentesLCH-noOUTLIERS.csv",
		  outputDir="data-ADDITIVE-noOUTLIERS")
#, nGenos=1000, nPheno=3)
gs_multi ("config-multi-ADDITIVE-noOUTLIERS.yml", "outputs-ADDITIVE-noOUTLIERS")
