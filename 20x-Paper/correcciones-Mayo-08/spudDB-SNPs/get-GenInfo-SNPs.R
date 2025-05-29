# Load required libraries
library(readxl)   # For reading .xls files
library(dplyr)    # For data manipulation
library(readr)    # For efficient CSV reading

# Get file paths from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript merge_snp_tables.R <input.xls> <input.csv>", call. = FALSE)
}

xls_file <- args[1]  # First argument: .xls file
csv_file <- args[2]  # Second argument: .csv file

# Read tables
table1 <- read_excel(xls_file)  # Read .xls (supports .xlsx too)
table2 <- read_csv(csv_file)    # Read .csv

# Merge tables by SNP column (inner join: keep only matching SNPs)
merged_table <- inner_join(
  table1 %>% select(SNP),       # Keep only SNP column from Table 1
  table2,                       # Keep all columns from Table 2
  by = "SNP"                    # Merge key
)

# Save merged table to a new CSV
output_file <- "merged_snps.csv"
write_csv(merged_table, output_file)

# Print summary
cat(sprintf(
  "Merged table saved to '%s'.\nInput SNPs: %d | Matched SNPs: %d\n",
  output_file,
  nrow(table1),
  nrow(merged_table)
))
