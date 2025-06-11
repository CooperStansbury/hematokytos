#!/usr/bin/env Rscript

# Utility: Color & Timestamp logging
log_msg <- function(msg, level = "INFO") {
  time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  color <- switch(level,
    "INFO"  = "\033[1;34m", # Blue
    "OK"    = "\033[1;32m", # Green
    "WARN"  = "\033[1;33m", # Yellow
    "ERROR" = "\033[1;31m", # Red
    "\033[0m"
  )
  reset <- "\033[0m"
  cat(sprintf("%s%s [%s] %s%s\n", color, time, level, msg, reset))
}

# Required Libraries
if (!requireNamespace("sceasy", quietly = TRUE)) {
  log_msg("Please install 'sceasy'", "ERROR")
  quit(status = 1)
}
if (!requireNamespace("reticulate", quietly = TRUE)) {
  log_msg("Please install 'reticulate'", "ERROR")
  quit(status = 1)
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  log_msg("Please install 'Seurat'", "ERROR")
  quit(status = 1)
}

library(sceasy)
library(reticulate)
library(Seurat)

# Set up Python environment for sceasy (customize if needed)
reticulate::use_condaenv('sceasy', required = TRUE)

# Argument Handling
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  log_msg("Usage: Rscript convert_seurat_to_anndata.R <input_rds_file> <output_h5ad_file>", "ERROR")
  quit(status = 1)
}

input_file <- args[1]
output_file <- args[2]

# Validate Input
if (!file.exists(input_file)) {
  log_msg(paste("Input file does not exist:", input_file), "ERROR")
  quit(status = 1)
}

# Ensure Output Directory Exists
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  log_msg(paste("Created output directory:", output_dir), "OK")
}

# Convert File
log_msg(sprintf("Converting '%s' to '%s'", basename(input_file), basename(output_file)), "INFO")
tryCatch({
  obj <- readRDS(input_file)
  # Check object type
  if (!("Seurat" %in% class(obj))) stop("Not a Seurat object")
  sceasy::convertFormat(obj, from = "seurat", to = "anndata", outFile = output_file)
  log_msg(sprintf("Saved: %s", output_file), "OK")
}, error = function(e) {
  log_msg(sprintf("Failed to convert %s: %s", basename(input_file), e$message), "ERROR")
  quit(status = 1)
})

log_msg("Conversion complete.", "INFO")
