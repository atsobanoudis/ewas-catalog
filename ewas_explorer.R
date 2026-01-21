# EWAS Catalog CpG Extractor
# Input: cpgs.txt
# Output: cpg_associations_full.json + summary CSV

if (!require("devtools", quietly = TRUE)) install.packages("devtools")
if (!require("ewascatalog", quietly = TRUE)) devtools::install_github("MRCIEU/ewascatalog-r")
if (!require("jsonlite", quietly = TRUE)) install.packages("jsonlite")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ewascatalog)
library(jsonlite)
library(dplyr)

# Read CpG list from file
cpg_file <- "cpgs.txt"
if (!file.exists(cpg_file)) {
  stop("cpgs.txt not found. Create a file with one CpG ID per line.")
}
cpg_list <- readLines(cpg_file) %>% trimws() %>% .[. != ""]

cat(sprintf("Loaded %d CpGs from %s\n", length(cpg_list), cpg_file))

# Query function
query_cpg <- function(cpg_id, p_threshold = 0.001) {
  tryCatch({
    result <- ewascatalog(cpg_id, "cpg")
    
    if (is.null(result) || nrow(result) == 0) {
      return(list(
        cpg_id = cpg_id,
        total_associations = 0,
        associations = list()
      ))
    }
    
    # Filter by p-value if column exists
    if ("p" %in% names(result) && is.numeric(result$p)) {
      result <- result %>% filter(p < p_threshold)
    }
    
    # Convert to clean list format
    associations <- result %>%
      mutate(across(where(is.factor), as.character)) %>%
      as.list() %>%
      lapply(function(x) if(length(x) == 1) x else list(x))
    
    list(
      cpg_id = cpg_id,
      total_associations = nrow(result),
      associations = result
    )
    
  }, error = function(e) {
    warning(sprintf("Error querying %s: %s", cpg_id, e$message))
    list(
      cpg_id = cpg_id,
      total_associations = 0,
      associations = list(),
      error = e$message
    )
  })
}

# Process all CpGs
cat("\nQuerying EWAS Catalog...\n")
results <- list()
pb <- txtProgressBar(min = 0, max = length(cpg_list), style = 3)

for (i in seq_along(cpg_list)) {
  cpg <- cpg_list[i]
  results[[cpg]] <- query_cpg(cpg)
  setTxtProgressBar(pb, i)
  Sys.sleep(0.3)  # Rate limiting
}
close(pb)

# Build metadata
metadata <- list(
  analysis_date = Sys.time(),
  total_cpgs = length(cpg_list),
  cpgs_with_results = sum(sapply(results, function(x) x$total_associations > 0)),
  database = "EWAS Catalog",
  p_threshold = 0.001
)

# Compile final output
output <- list(
  metadata = metadata,
  results = results
)

# Save JSON
json_file <- "cpg_associations_full.json"
write_json(output, json_file, pretty = TRUE, auto_unbox = TRUE)
cat(sprintf("\nSaved: %s\n", json_file))

# Create summary CSV
summary_df <- data.frame(
  cpg_id = names(results),
  total_associations = sapply(results, function(x) x$total_associations),
  has_data = sapply(results, function(x) x$total_associations > 0)
)

csv_file <- "cpg_associations_summary.csv"
write.csv(summary_df, csv_file, row.names = FALSE)
cat(sprintf("Saved: %s\n", csv_file))

# Print summary
cat(sprintf("\n=== SUMMARY ===\n"))
cat(sprintf("Total CpGs queried: %d\n", metadata$total_cpgs))
cat(sprintf("CpGs with associations: %d\n", metadata$cpgs_with_results))
cat(sprintf("CpGs with no data: %d\n", sum(!summary_df$has_data)))
cat(sprintf("Total associations: %d\n", sum(summary_df$total_associations)))