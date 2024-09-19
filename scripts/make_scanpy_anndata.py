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
    annotations_path = sys.argv[2]
    output_path = sys.argv[3]

    cell_type_path = sys.argv[2]
    output_path = sys.argv[3]
    
    # Load the data
    print(f"\nLoading data from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    
    # Load cell types
    print(f"\nLoading annotations from {annotations_path}...")
    df = load_cell_types(annotations_path)
    
    print(f"\nMerging annotations...")
    adata.obs = pd.merge(
        adata.obs,
        df,
        how='left',
        left_index=True,
        right_index=True,
    )
    
    
    
    
    
    
    
#     # Extract the observation metadata
#     df = adata.obs.copy()
#     df = df.reset_index(drop=False, names='obs_index')

#     # Save the data
#     print(f"\nSaving metadata to {output_path}...")
#     df.to_csv(output_path, index=False)
#     print("Done!")
    
    

 
    
    
    
    