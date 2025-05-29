
source ("lglib14.R")
#-------------------------------------------------------------
#-------------------------------------------------------------
main <- function () {
	args     = commandArgs (trailingOnly=T)
	fenoFile = "fenotiposColor-ComponentesLCH.csv"
    genoFile = "genotipo-AndigenaCCC-ClusterCall2020-MATRIX.csv"

    simplify_tables (fenoFile, genoFile, 100, 10)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
simplify_tables <- function(traits_file, markers_file, N, M) {
  # Read both CSV files
  traits <- read.csv(traits_file, stringsAsFactors = FALSE, check.names=F)
  markers <- read.csv(markers_file, stringsAsFactors = FALSE, check.names=F)

  # Order reg numbers alphabetically from first column of traits
  ordered_regs <- sort(traits[[1]])
  selected_regs <- head(ordered_regs, N)

  # Subset traits table
  trait_subset <- traits[traits[[1]] %in% selected_regs, 1:4]  # First column + 3 trait columns

  # Subset markers table
  markersCols = colnames(markers)
  view (markersCols)
  view (selected_regs)


  marker_subset <- markers[1:M, c(TRUE, colnames(markers)[-1] %in% selected_regs)]

  # Write output
  write.csv(trait_subset, sub(".csv", "_TEST.csv", traits_file), row.names = FALSE)
  write.csv(marker_subset, sub(".csv", "_TEST.csv", markers_file), row.names = FALSE)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
main ()
