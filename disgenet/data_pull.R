# Load required packages, installing if necessary.
if (!requireNamespace("disgenet2r", quietly = TRUE)) {
  install.packages("disgenet2r")
}
if (!requireNamespace("dotenv", quietly = TRUE)) {
  install.packages("dotenv")
}
library(disgenet2r)
library(dotenv)

# Load environment variables
load_dot_env()

# --- Configuration ---
# Define the path to the file containing the list of gene symbols.
gene_list_file <- "data/gene_list.txt"

# Read the list of genes from the text file.
if (file.exists(gene_list_file)) {
  genes_of_interest <- readLines(gene_list_file)
  # Remove any empty lines or whitespace around gene symbols
  genes_of_interest <- trimws(genes_of_interest)
  genes_of_interest <- genes_of_interest[genes_of_interest != ""]
} else {
  stop("Error: Gene list file not found at ", gene_list_file)
}

# Define output file paths
gda_output_file <- "data/gene_disease_associations.rds"
gea_output_file <- "data/gene_evidence_associations.rds"


# --- Data Fetching ---
# Check if the data has already been fetched. If not, call the API.

# 1. Gene-Disease Associations
if (!file.exists(gda_output_file)) {
  cat("Fetching gene-disease associations for:",
      paste(genes_of_interest, collapse = ", "), "\n")
  gda_results <- gene2disease(gene = genes_of_interest, database = "CURATED")
  cat("Saving gene-disease results to", gda_output_file, "\n")
  saveRDS(gda_results, file = gda_output_file)
} else {
  cat("Gene-disease association data already exists:", gda_output_file, "\n")
}

# 2. Gene-Evidence Associations
if (!file.exists(gea_output_file)) {
  cat("Fetching gene-evidence associations for:",
      paste(genes_of_interest, collapse = ", "), "\n")
  gea_results <- gene2evidence(gene = genes_of_interest, database = "CURATED")
  cat("Saving gene-evidence results to", gea_output_file, "\n")
  saveRDS(gea_results, file = gea_output_file)
} else {
  cat("Gene-evidence association data already exists:", gea_output_file, "\n")
}

cat("Data pull complete.\n")