import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as an
from scipy.sparse import csr_matrix

# local imports
import utils as ut

sc.settings.verbosity = 3  


def skeletonize(adata):
    """
    Reduces the size of an AnnData object by removing non-essential data.

    Args:
        adata: An AnnData object.

    Returns:
        A streamlined AnnData object.
    """

    # Remove unneeded annotations
    for attr in ['uns', 'obsm', 'obsp', 'varm']:
        if hasattr(adata, attr):
            delattr(adata, attr)

    # Handle layers
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
    elif 'raw_counts' in adata.layers:
        adata.X = adata.layers['raw_counts'].copy()
    del adata.layers

    # Make sparse
    adata.X = csr_matrix(adata.X)

    # Handle var
    if 'gene_name' in adata.var.columns:
        adata.var_names = adata.var['gene_name']
    del adata.var

    # Handle obs
    del adata.obs

    return adata



        
if __name__ == "__main__":
    adata_path = sys.argv[1]
    min_cells = int(sys.argv[2])
    min_genes = int(sys.argv[3])
    target_sum = int(eval(sys.argv[4]))
    token_table = sys.argv[5]
    output_path = sys.argv[6]

    # Load the data
    print(f"\nLoading data from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    print("Memory usage after loading:")
    sc.logging.print_memory_usage()

    # Skeletonize
    print("\nSkeletonizing data...")
    adata = skeletonize(adata)
    print("Memory usage after skeletonizing:")
    sc.logging.print_memory_usage()
    
    print("\nFiltering genes by GeneFormer corpus...")
    # remove genes not in the corpus
    tokens = pd.read_csv(token_table)
    adata = ut.filter_genes_and_tokens(adata, tokens)
    
    # Filtration
    print("\nFiltering cells and genes...")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    adata.layers["counts"] = adata.X.copy()

    # Normalization
    print("\nNormalizing and transforming data...")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # Save the data
    print(f"\nSaving processed data to {output_path}...")
    adata.write(output_path)
    print("Done!")
    
    

 
    
    
    
    