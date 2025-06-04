# Pipelines

This directory houses the workflow scripts used throughout the project. Each subfolder contains the configuration and helper scripts for a specific analysis pipeline.

- `cc_pipeline/` – Nextflow workflow for processing cell cycle datasets
- `bm_pipeline/` – Nextflow workflow for bone marrow single-cell data
- `hsc_pipeline/` – Nextflow workflow for direct reprogramming experiments
- `velocyto_pipeline/` – Snakemake pipeline to produce Velocyto spliced/unspliced counts
- `DeepCycle/` – Scripts to run the DeepCycle model
- `CABYBARA/` – Capybara-based cell identity assignment scripts
- `scMINER/` – Gene regulatory network inference using SJARACNe

Each pipeline folder includes a README describing how to run the workflow.
