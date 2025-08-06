# Resource Files

This directory collects gene sets, metadata, and auxiliary files used
across pipelines and notebooks.

## Gene sets
- `CellMarker_2024.txt` – Canonical markers from the 2024 CellMarker release
- `CellMarker_Augmented_2021.txt` – Extended CellMarker table from 2021
- `GO_Biological_Process_2023.txt` – Gene Ontology gene sets
- `PanglaoDB_Augmented_2021.txt` – PanglaoDB marker genes
- `Tabula_Sapiens.txt` – Cell type markers extracted from Tabula Sapiens
- `allTFs_hg38.txt` – List of human transcription factors
- `human_cell_cycle_genes.csv` – Curated cell-cycle genes
- `kamimoto_genes.csv` – Marker genes reported by Kamimoto et al.
- `literature_genes.csv` – Additional genes collated from the literature
- `quiescence_markers.csv` – Gene set used for quiescence scoring
- `regev_lab_cell_cycle_genes.txt` – Regev lab cell-cycle genes
- `seurat_genes.csv` – Marker genes from the Seurat package

## Metadata and configuration
- `cell_labels.yaml` – Cell ontology annotations and groupings
- `cell_map.yaml` – Mapping between ontology labels and common names
- `cell_type_annotation.csv` – Example cell type annotations
- `clean_marker_genes.csv` – Filtered marker gene table
- `gene_names.tsv.gz` – Alias table for gene symbols
- `geneformer_params.yaml` – Configuration for the Geneformer model
- `hsc_cell_types.yaml` – Groupings for HSC-related cell types
- `cl.owl` – Cell Ontology in OWL format

## Utility files
- `make_marker_gene_df.ipynb` – Notebook to combine marker lists
- `example_marker_data.csv` – Sample dataset for documentation examples
