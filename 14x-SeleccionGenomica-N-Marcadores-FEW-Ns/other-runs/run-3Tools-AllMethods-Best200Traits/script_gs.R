#!/usr/bin/env Rscript
#detach("package:ppGS", unload=TRUE)
library (ppGS)

gs_multi ('params-nMarkers.yml', "outputs")
