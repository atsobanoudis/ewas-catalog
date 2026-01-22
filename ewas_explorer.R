# EWAS Atlas CpG Extractor

if (!require("httr", quietly = TRUE)) install.packages("httr")
if (!require("jsonlite", quietly = TRUE)) install.packages("jsonlite")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")

library(httr)
library(jsonlite)
library(dplyr)

API_ROOT <- "https://ngdc.cncb.ac.cn/ewas/rest"

# Extract probe data
extract_probe <- function(probe_id) {
  url <- sprintf("%s/probe?probeId=%s", API_ROOT, probe_id)
  
  tryCatch({
    response <- GET(url, timeout(30))
    
    if (status_code(response) != 200) {
      return(list(probe_id = probe_id, error = sprintf("HTTP %d", status_code(response))))
    }
    
    result <- fromJSON(content(response, "text", encoding = "UTF-8"))
    
    if (result$code != 0 || is.null(result$data)) {
      return(list(probe_id = probe_id, error = "No data returned"))
    }
    
    d <- result$data
    
    # Extract correlation counts
    assoc <- d$associationList
    corr_counts <- if (!is.null(assoc) && nrow(assoc) > 0) {
      table(assoc$correlation)
    } else {
      c(hyper = 0, hypo = 0, NR = 0)
    }
    
    # Extract unique genes
    genes <- if (!is.null(d$relatedTranscription) && nrow(d$relatedTranscription) > 0) {
      unique(d$relatedTranscription$geneName)
    } else {
      character(0)
    }
    
    # Extract traits
    traits <- if (!is.null(assoc) && nrow(assoc) > 0) {
      unique(assoc$trait)
    } else {
      character(0)
    }
    
    # Build result
    list(
      probe_id = probe_id,
      chr = d$chrHg19,
      pos = d$posHg19,
      cpg_island = d$cpgIsland,
      n_studies = if (!is.null(assoc)) nrow(assoc) else 0,
      correlations = list(
        hyper = as.integer(corr_counts["pos"] %||% 0),
        hypo = as.integer(corr_counts["neg"] %||% 0),
        NR = as.integer(corr_counts["NR"] %||% 0)
      ),
      related_genes = genes,
      related_traits = traits,
      studies = if (!is.null(assoc)) assoc else data.frame()
    )
    
  }, error = function(e) {
    list(probe_id = probe_id, error = e$message)
  })
}

# Process multiple CpGs
process_cpgs <- function(cpg_list) {
  results <- list()
  pb <- txtProgressBar(min = 0, max = length(cpg_list), style = 3)
  
  for (i in seq_along(cpg_list)) {
    cpg <- cpg_list[i]
    results[[cpg]] <- extract_probe(cpg)
    setTxtProgressBar(pb, i)
    Sys.sleep(0.1)  # Polite rate limiting
  }
  
  close(pb)
  return(results)
}

# Main execution
cat("EWAS Atlas CpG Extractor\n\n")

cpg_file <- "cpgs.txt"
if (!file.exists(cpg_file)) {
  stop("cpgs.txt not found")
}

cpg_list <- readLines(cpg_file) %>% trimws() %>% .[. != ""]
cat(sprintf("Loaded %d CpGs\n\n", length(cpg_list)))

cat("Querying EWAS Atlas API...\n")
results <- process_cpgs(cpg_list)

# Build output
output <- list(
  metadata = list(
    analysis_date = Sys.time(),
    total_cpgs = length(cpg_list),
    api_source = "EWAS Atlas",
    api_root = API_ROOT
  ),
  results = results
)

# Save JSON
write_json(output, "cpg_associations_full.json", pretty = TRUE, auto_unbox = TRUE)
cat("\nSaved: cpg_associations_full.json\n")

# Create summary CSV
summary_df <- lapply(results, function(r) {
  data.frame(
    cpg_id = r$probe_id,
    chr = r$chr %||% NA,
    pos = r$pos %||% NA,
    cpg_island = r$cpg_island %||% NA,
    n_studies = r$n_studies %||% 0,
    hyper = r$correlations$hyper %||% 0,
    hypo = r$correlations$hypo %||% 0,
    NR = r$correlations$NR %||% 0,
    n_genes = length(r$related_genes %||% character(0)),
    genes = paste(r$related_genes %||% character(0), collapse = ", "),
    n_traits = length(r$related_traits %||% character(0)),
    error = !is.null(r$error)
  )
}) %>% bind_rows()

write.csv(summary_df, "cpg_associations_summary.csv", row.names = FALSE)
cat("Saved: cpg_associations_summary.csv\n")

cat(sprintf("\nProcessed %d CpGs\n", nrow(summary_df)))
cat(sprintf("With associations: %d\n", sum(summary_df$n_studies > 0)))
cat(sprintf("Errors: %d\n", sum(summary_df$error)))