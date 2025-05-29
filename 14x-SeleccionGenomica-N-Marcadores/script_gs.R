#!/usr/bin/env Rscript

Sys.setenv(RENV_PROJECT = "/home/lg/BIO/agrosavia/genetica-color-papa/14x-SeleccionGenomica-N-Marcadores")
source(file.path(Sys.getenv("RENV_PROJECT"), "renv/activate.R"))
#library (ppGS)
devtools::load_all(file.path(Sys.getenv("RENV_PROJECT"), "lib/package-ppGS"))

gs_multi ('params-nMarkers.yml', "outputs")
