import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as an
from scipy.sparse import csr_matrix

sc.settings.verbosity = 3  


if __name__ == "__main__":
    adata_path = sys.argv[1]
    obs_path = sys.argv[2]
    var_path = sys.argv[3]

    # Load the data
    print(f"\nLoading data from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    
    # extract QC metrics
    obs, var = sc.pp.calculate_qc_metrics(
        adata,
        layer='counts',
        inplace=False
    )

    obs = obs.reset_index(drop=False, names='obs_index')
    var = var.reset_index(drop=False, names='var_index')

    print(f"{obs.shape=}")
    print(f"{var.shape=}")
    
    # Save the data
    print(f"\nSaving obs metrics to {obs_path}...")
    obs.to_csv(obs_path, index=False)
    
    print(f"\nSaving var metrics to {var_path}...")
    var.to_csv(var_path, index=False)
    
    print("Done!")
    
    

 
    
    
    
    