#!/usr/bin/env Rscript

# Disable package startup messages
suppressPackageStartupMessages({
  library(AnnotationHub)
})

# Nice console prints
message("▶ Initializing AnnotationHub...")
hub <- AnnotationHub()

message("▶ Querying alternative splicing events for hg38...")
events <- suppressMessages(
  query(hub, "alternativeSplicingEvents.hg38_V2")
)[[1]]

outdir <- "/nfs/turbo/umms-indikar/shared/projects/HSC/data/alternative_splicing"
message("▶ Ensuring output directory exists: ", outdir)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  message("   • Created directory.")
} else {
  message("   • Directory already exists.")
}

n_events <- length(events)
message("▶ Found ", n_events, " event type", if (n_events != 1) "s" else "", ".")

# Loop and write
for (nm in names(events)) {
  message("\n→ Processing event type: ", nm)
  df <- events[[nm]]
  
  # Detect list‐columns
  list_cols <- vapply(df, is.list, logical(1))
  if (any(list_cols)) {
    message("   • Flattening ", sum(list_cols), " list-columns into strings.")
    df[list_cols] <- lapply(df[list_cols], function(col) {
      vapply(col, function(x) paste(x, collapse = ";"), character(1))
    })
  }
  
  # Sanitize name and write CSV
  fname <- paste0(gsub("[^A-Za-z0-9_]", "_", nm), ".csv")
  fpath <- file.path(outdir, fname)
  message("   • Writing CSV: ", fpath)
  write.csv(df, file = fpath, row.names = FALSE)
}

message("\n✔ All done. CSV files are in: ", outdir)
