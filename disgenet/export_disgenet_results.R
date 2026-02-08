# Load required packages
if (!requireNamespace("disgenet2r", quietly = TRUE)) {
  stop("disgenet2r package is required but not installed.")
}
library(disgenet2r)

# Function to safely extract qresult and flatten list-columns
extract_results <- function(rds_file) {
  if (!file.exists(rds_file)) {
    warning(paste("File not found:", rds_file))
    return(NULL)
  }
  obj <- readRDS(rds_file)
  if (.hasSlot(obj, "qresult")) {
    df <- obj@qresult
    
    # Flatten list-columns to strings to prevent write.csv errors
    list_cols <- sapply(df, is.list)
    if (any(list_cols)) {
      for (col_name in names(df)[list_cols]) {
        df[[col_name]] <- sapply(df[[col_name]], function(x) {
          if (is.null(x) || length(x) == 0) return("")
          paste(as.character(x), collapse = " | ")
        })
      }
    }
    return(df)
  } else {
    warning(paste("No 'qresult' slot found in", rds_file))
    return(NULL)
  }
}

# --- 1. Load Input Data ---
gene_list_file <- "data/gene_list.txt"
if (file.exists(gene_list_file)) {
  input_genes <- readLines(gene_list_file)
  input_genes <- trimws(input_genes)
  input_genes <- input_genes[input_genes != ""]
  # Handle potential duplicates in input
  input_genes <- unique(input_genes)
  cat("Loaded", length(input_genes), "genes from input list.
")
} else {
  stop("Gene list file not found.")
}

# --- 2. Process Gene-Disease Associations (GDA) ---
gda_file <- "data/gene_disease_associations.rds"
cat("
Processing Gene-Disease Associations...
")
gda_df <- extract_results(gda_file)

if (!is.null(gda_df)) {
  # Export to CSV
  gda_csv_file <- "data/disgenet_gda.csv"
  write.csv(gda_df, gda_csv_file, row.names = FALSE)
  cat("  Saved detailed GDA to:", gda_csv_file, "
")
  
  # Check for missing genes
  # Note: The column name for gene symbol in disgenet2r results is typically 'gene_symbol'
  if ("gene_symbol" %in% colnames(gda_df)) {
    found_genes_gda <- unique(gda_df$gene_symbol)
    missing_genes_gda <- setdiff(input_genes, found_genes_gda)
    
    cat("  Genes with associations:", length(found_genes_gda), "
")
    cat("  Genes without associations:", length(missing_genes_gda), "
")
  } else {
    warning("  'gene_symbol' column not found in GDA results. Cannot compute missing genes.")
  }
}

# --- 3. Process Gene-Evidence Associations (GEA) ---
gea_file <- "data/gene_evidence_associations.rds"
cat("
Processing Gene-Evidence Associations...
")
gea_df <- extract_results(gea_file)

if (!is.null(gea_df)) {
  # Export to CSV
  gea_csv_file <- "data/disgenet_gea.csv"
  write.csv(gea_df, gea_csv_file, row.names = FALSE)
  cat("  Saved detailed GEA to:", gea_csv_file, "
")
}

# --- 4. Analyze Missing Genes ---
# We focus on the GDA results for the "missing" analysis as that's the primary association map.
if (!is.null(gda_df) && exists("missing_genes_gda")) {
  
  # Create a summary dataframe for missing genes
  # Since we can't easily distinguish "Not in DB" vs "No Associations" without querying the DB for every gene individually
  # or parsing the API logs, we will list them as "No Association Found (in CURATED DB)".
  
  missing_df <- data.frame(
    gene_symbol = missing_genes_gda,
    status = "No association found in CURATED database"
  )
  
  missing_csv_file <- "data/genes_missing_associations.csv"
  write.csv(missing_df, missing_csv_file, row.names = FALSE)
  cat("
Saved list of genes with no associations to:", missing_csv_file, "
")
  
  # Print a preview of missing genes
  if (length(missing_genes_gda) > 0) {
    cat("
Example genes with no associations:
")
    print(head(missing_genes_gda, 10))
  }
}

cat("
Processing complete.
")
