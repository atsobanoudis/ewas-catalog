# Load required packages
if (!requireNamespace(c("disgenet2r", "devtools"), quietly = TRUE)) {
  install.packages("devtools")
  library(devtools)
  install_gitlab("medbio/disgenet2r")
}
library(disgenet2r)

# --- Configuration ---
# Define the input file path. This should match the output from data_pull.R
gda_input_file <- "data/gene_disease_associations.rds"

# --- Analysis ---
# Check if the data file exists before trying to load it.
if (file.exists(gda_input_file)) {
  # Load the saved gene-disease association results
  cat("Loading gene-disease data from", gda_input_file, "\n")
  gda_results <- readRDS(gda_input_file)

  # --- Visualization ---
  # Create a heatmap of the gene-disease associations by disease class.
  # The plot function is a built-in method for DataGeNET.DGN objects.
  cat("Generating heatmap...\n")

  # The plot function may open an interactive window or save to a file
  # depending on the R environment.
  plot(
    gda_results,
    type = "Heatmap",
    class = "DiseaseClass",
    nchars = 60,
    interactive = TRUE
  )

  cat("Analysis script finished.\n")

} else {
  cat("Error: Data file not found at", gda_input_file, "\n")
  cat("Please run data_pull.R first to download the data.\n")
}