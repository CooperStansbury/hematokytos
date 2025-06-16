# Analysis Notebooks

This folder contains the bulk of the exploratory analyses performed in the project.
Notebooks are organised by topic with additional subdirectories for larger groups
of related notebooks.

## Subdirectories
- `capybara/` – Experiments using the Capybara label transfer tool
- `cell_assign/` – Notebooks applying the CellAssign algorithm to our data
- `cell_cycle/` – Evaluation of cell-cycle state including DeepCycle results
- `geneformer/` – Exploration of the Geneformer language model
- `isoforms/` – Isoform-level expression analysis
- `reference_analysis/` – Building and interrogating the reference atlas
- `tf_networks/` – Transcription factor network inference with scMINER and others
- `trajectory_inference/` – Pseudotime and RNA velocity analyses

## Selected notebooks
- `DEG_our_data.ipynb` – Differential expression overview
- `anndata_summary.ipynb` – Inspect AnnData objects and metadata
- `distance_analysis.ipynb` – Compare sample distances
- `marker_gene_expression.ipynb` – Plot canonical marker genes
- `scMINER_prepare.ipynb` – Formatting data for network inference
- `sequencing_summary.ipynb` – Summaries of sequencing metrics
