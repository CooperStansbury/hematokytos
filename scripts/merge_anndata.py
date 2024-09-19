import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as an
from scipy.sparse import csr_matrix

# local import
import utils as ut

sc.settings.verbosity = 3  


if __name__ == "__main__":
    output_path = sys.argv[1]
    token_table = sys.argv[2]
    file_list = sys.argv[3:]
    
    adata_list = {}
    
    for adata_path in file_list:
        basename = os.path.basename(adata_path)
        sample_id = basename.replace(".h5ad", "")
        
        print(f"\nLoading data from {adata_path}...")
        adata = sc.read_h5ad(adata_path)
        print("Memory usage after loading:")
        sc.logging.print_memory_usage()
        
        adata.obs['dataset'] = sample_id
        
        adata_list[sample_id] = adata

    print(f"\nMerging all AnnData objects...")
    
    adata = an.concat(
        adata_list, 
        axis=0, # align on obs
        label='dataset', 
        index_unique="_", 
        join="outer",
        merge='same',
    )

    print("Memory usage after merging:")
    sc.logging.print_memory_usage()
    
    print("Standard processing:")
    adata = ut.quick_process_anndata(adata)
    
    print("\nFiltering genes by GeneFormer corpus...")
    # remove genes not in the corpus
    tokens = pd.read_csv(token_table)
    adata = ut.filter_genes_and_tokens(adata, tokens)
    
    print("Memory usage after processing:")
    sc.logging.print_memory_usage()
    
    print("\nAnnData:")
    print(adata)

    # Save the data
    print(f"\nSaving processed data to {output_path}...")
    adata.write(output_path)
    print("Done!")
    
    
    
    