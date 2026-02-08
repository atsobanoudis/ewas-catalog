# Load required packages
if (!requireNamespace("disgenet2r", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  library(devtools)
  install_gitlab("medbio/disgenet2r")
}
library(disgenet2r)

# --- Configuration ---
# Define the input file path
gda_input_file <- "data/gene_disease_associations.rds"
output_file <- "viz/heatmap_disgenet2r.png"

# Ensure output directory exists
if (!dir.exists("viz")) {
  dir.create("viz")
}

# --- Analysis ---
if (file.exists(gda_input_file)) {
  cat("Loading gene-disease data from", gda_input_file, "
")
  gda_results <- readRDS(gda_input_file)

  # --- Visualization ---
  cat("Generating static heatmap (300 DPI)...
")
  
  # Open PNG graphics device
  png(filename = output_file, width = 12, height = 10, units = "in", res = 300)
  
  # Generate Plot using disgenet2r
  # Note: The 'plot' function in disgenet2r typically returns a ggplot object or prints to device.
  # We wrap it in tryCatch to handle potential internal plotting issues.
  tryCatch({
      plot(
        gda_results,
        type = "Heatmap",
        class = "DiseaseClass",
        nchars = 60,
        interactive = FALSE
      )
  }, error = function(e) {
      cat("Error plotting heatmap:", e$message, "
")
  })
  
  # Close device
  dev.off()

  cat("Heatmap saved to:", output_file, "
")

} else {
  cat("Error: Data file not found at", gda_input_file, "
")
}
