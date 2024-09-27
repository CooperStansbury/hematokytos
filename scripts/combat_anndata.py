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


if __name__ == "__main__":
    adata_path = sys.argv[1]
    output_path = sys.argv[2]
    
    print(f"\nLoading data from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    print("Memory usage after loading:")
    sc.logging.print_memory_usage()
      
    print(f"\nRunning batch correction using combat...")

    # combat batch correction
    sc.pp.combat(
        adata, 
        key='dataset',
        covariates=['total_counts', 'n_genes'],
    )
    
    # make matrix sparse, reduce the precision
    adata.X = csr_matrix(adata.X.astype('float32'))
    
    print("Memory usage after combat:")
    sc.logging.print_memory_usage()
    
    print("Standard processing:")
    adata = ut.quick_process_anndata(adata)
    
    print("Memory usage after processing:")
    sc.logging.print_memory_usage()
    
    print("\nAnnData:")
    print(adata)

    # Save the data
    print(f"\nSaving processed data to {output_path}...")
    adata.write(output_path)
    print("Done!")
    