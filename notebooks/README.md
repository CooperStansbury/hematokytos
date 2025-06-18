# Analysis Notebooks

This directory contains the bulk of our exploratory analyses and figure preparation workflows.
Notebooks are grouped by topic with subfolders collecting related experiments. Each subfolder contains a small README describing its contents.

## Subdirectories
- `benchmarks/` – Scripts comparing datasets and integration strategies
- `bone_marrow/` – Quality control for bone marrow samples
- `capybara/` – Experiments using the Capybara label transfer tool
- `cell_assign/` – Example runs of the CellAssign algorithm
- `cell_cycle/` – Estimating cell-cycle state and related QC
- `geneformer/` – Exploration of the Geneformer language model
- `isoforms/` – Isoform-level expression analysis
- `reference_analysis/` – Building and interrogating the reference atlas
- `summary/` – General summaries and dataset overviews
- `tf_networks/` – Transcription factor network inference with scMINER and others
- `trajectory_inference/` – Pseudotime and RNA velocity analyses

## Selected notebooks
- `DEG_our_data.ipynb` – Differential expression overview
- `anndata_summary.ipynb` – Inspect AnnData objects and metadata
- `distance_analysis.ipynb` – Compare sample distances
- `marker_gene_expression.ipynb` – Plot canonical marker genes
- `scMINER_prepare.ipynb` – Formatting data for network inference
- `sequencing_summary.ipynb` – Summaries of sequencing metrics
- `nearest_neighbors.ipynb` – Quick nearest neighbor visualisations
- `make_colorbar.ipynb` – Utility for consistent plotting colours
