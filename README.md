# Hematokytos

Hematokytos ("blood cell") collects pipelines and utilities for gathering,
integrating, and analysing single-cell data from direct reprogramming
experiments alongside native hematopoietic cell types. Workflows cover
preprocessing, quality control, RNA velocity, and network analysis with a
focus on hematopoietic stem and progenitor cells.

## Repository Structure

- `pipelines/` – Collection of workflow scripts
  - `cc_pipeline/` – Nextflow pipeline for cell cycle samples
  - `bm_pipeline/` – Nextflow pipeline for bone marrow data
  - `hsc_pipeline/` – Nextflow pipeline for direct reprogramming experiments
  - `velocyto_pipeline/` – Snakemake workflow for running Velocyto
  - `DeepCycle/` – Scripts to run the DeepCycle model
  - `CABYBARA/` – Capybara cell identity assignment scripts
  - `scMINER/` – Network inference with SJARACNe
- `notebooks/` – Jupyter notebooks for exploratory analysis and figure generation
  (see `notebooks/README.md` for an overview)
- `reference_atlas/` – Scripts for collecting and annotating reference datasets
- `resources/` – Gene sets, metadata and supporting files
- `NextCell/` – Exploration of the NextCell algorithm
- `results/` – Example output tables
- `utils/` – Small helper modules

Each pipeline directory contains its own README with usage instructions.

