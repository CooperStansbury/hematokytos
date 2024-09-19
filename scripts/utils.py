import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import yaml
import scanpy.external as sce
import anndata as an
from scipy.sparse import csr_matrix

sc.settings.verbosity = 3 


def load_yaml(file_path):
    """Loads a YAML file and returns its contents as a Python dictionary.

    Args:
    file_path: The path to the YAML file.

    Returns:
    A Python dictionary representing the YAML data.
    """
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
    return data


def load_annotations(fpath):
    """
    Loads annotations data from a CSV file, processes it, and returns the resulting DataFrame.

    Args:
        fpath (str): The file path to the CSV file containing annotations data.

    Returns:
        pd.DataFrame: The processed DataFrame with 'cell_id' as the index.
    """
    df = pd.read_csv(fpath)
    df['cell_id'] = df['obs_index'].astype(str) + "_" + df["dataset"]
    df = df.drop(columns='dataset')
    df = df.set_index('cell_id')
    return df


def filter_genes_and_tokens(adata, tokens):
    """Filters an AnnData object and a token DataFrame based on mutual gene presence.

    Args:
        adata: An AnnData object containing gene expression data.
        tokens: A DataFrame with 'gene_name' as a column and potentially other token-related information.

    Returns:
        The filtered AnnData object with token information merged into `.var`.
    """

    # Filter adata by genes present in tokens
    print("Filtering adata based on genes in tokens...")
    adata = adata[:, adata.var.index.isin(tokens['gene_name'].values)].copy()
    print(f"Retained {adata.shape[1]} genes in adata that are also present in tokens.")
    
    # Filter tokens by genes present in adata
    print("Filtering tokens based on genes in adata...")
    tokens = tokens[tokens['gene_name'].notna()]
    tokens = tokens.sort_values(by=['gene_name', 'gene_id'])
    tokens = tokens.drop_duplicates(subset='gene_name')
    tokens = tokens[tokens['gene_name'].isin(adata.var_names)]
    print(f"Retained {len(tokens)} tokens corresponding to genes in adata.")
    
    # Set 'gene_name' as the index in tokens for easier merging
    print("Setting 'gene_name' as index in tokens...")
    tokens = tokens.set_index('gene_name')

    # Combine adata.var and tokens, aligning on 'gene_name'
    print("Merging token information into adata.var...")
    adata.var = pd.concat([adata.var, tokens], ignore_index=False, axis=1)
    print("Merge completed!")

    return adata 



def quick_process_anndata(adata):
    """
    Performs standard, simple processing steps on an AnnData object.

    This function executes a series of common analysis steps on an AnnData object, 
    assuming that the data in `adata.X` is already normalized and log-transformed.

    Steps included:

    1. Calculate quality control (QC) metrics for each cell.

    Args:
        adata: An AnnData object containing normalized, log-transformed data in `adata.X`.

    Returns:
        The processed AnnData object with QC metrics
    """
    # QC metrics calculation
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    return adata

    

def process_anndata(adata):
    """
    Performs standard processing steps on an AnnData object.

    This function executes a series of common analysis steps on an AnnData object, 
    assuming that the data in `adata.X` is already normalized and log-transformed.

    Steps included:

    1. Calculate quality control (QC) metrics for each cell.
    2. Identify highly variable genes across batches (specified by 'dataset').
    3. Perform principal component analysis (PCA) on highly variable genes.
    4. Integrate data across batches using Scanorama.
    5. Compute a nearest neighbor graph.
    6. Calculate UMAP embeddings for visualization.
    7. Cluster cells using the Leiden algorithm.

    Args:
        adata: An AnnData object containing normalized, log-transformed data in `adata.X`.

    Returns:
        The processed AnnData object with QC metrics, highly variable genes, PCA 
        results, Scanorama integration, neighbor graph, UMAP embeddings, and Leiden 
        clusters.
    """

    # QC metrics calculation
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # Highly variable genes identification
    sc.pp.highly_variable_genes(adata, batch_key='dataset')

    # PCA
    sc.tl.pca(adata, mask_var='highly_variable')

    # Scanorama integration
    sce.pp.scanorama_integrate(adata, key='dataset')

    # Neighbor graph computation
    sc.pp.neighbors(adata)

    # UMAP embedding
    sc.tl.umap(adata, min_dist=0.5)

    # Leiden clustering
    sc.tl.leiden(adata, resolution=0.5)

    return adata
