# Pipelines

Workflow scripts used throughout the project live here. Each subdirectory
contains configuration files and helper scripts for a specific analysis
pipeline and includes its own README with usage details.

- `cc_pipeline/` – Nextflow workflow for processing cell cycle datasets
- `bm_pipeline/` – Nextflow workflow for bone marrow single-cell data
- `hsc_pipeline/` – Nextflow workflow for direct reprogramming experiments
- `velocyto_pipeline/` – Snakemake pipeline to produce Velocyto spliced/unspliced counts
- `DeepCycle/` – Scripts to run the DeepCycle model
- `CABYBARA/` – Capybara-based cell identity assignment scripts
- `scMINER/` – Gene regulatory network inference using SJARACNe
