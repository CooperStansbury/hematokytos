import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as an
from datasets import Dataset
from scipy.sparse import csr_matrix

sc.settings.verbosity = 3  


# local imports 
import utils as ut



def tokenize(adata, chunk_size=1000, max_tokens=2048, pad_token=0):
    """
    Tokenizes anndata objects for GeneFormer processing.

    Args:
        adata: Anndata object containing gene expression data.
        chunk_size: Number of cells to process at once (for memory efficiency).
        max_tokens: Maximum number of tokens to include per cell.

    Assumptions:
        1. adata.X values are normalized appropriately.
        2. Genes are already subset to the GeneFormer corpus with token IDs in adata.var.
        3. all needed annotations are merged into adata.obs 

    Returns:
        DataFrame with 'cell_id', 'input_ids', 'length', and 'total_length' columns.
    """

    token_ids = adata.var['token_id'].values
    result = []

    total_chunks = adata.shape[0] // chunk_size + 1  # Calculate total chunks for progress

    for i, start in enumerate(range(0, adata.shape[0], chunk_size)):
        end = start + chunk_size
        adata_chunk = adata[start:end, :]

        X = adata_chunk.to_df()
        X.columns = token_ids
    
        for cell_id, row in X.iterrows():
            ranks = row[row > 0].rank(method='first').astype(int).sort_values()
            input_ids = ranks.head(max_tokens).index.to_list()

            # Pad with pad_token if necessary
            padding_length = max_tokens - len(input_ids)
            if padding_length > 0:
                input_ids += [pad_token] * padding_length
            
            # get obs metadata
            new_row = adata_chunk.obs.loc[cell_id, ].to_dict()
            
            new_row['cell_id'] = cell_id
            new_row['input_ids'] = input_ids
            new_row['length'] = len(input_ids)
            new_row['total_length'] = len(ranks)
            result.append(new_row)
            
        # Print progress update
        if (i + 1) % 10 == 0 or i + 1 == total_chunks:  # Print every 10 chunks or on the last chunk
            print(f"Processed {min(end, adata.shape[0])} out of {adata.shape[0]} cells ({(i + 1) / total_chunks:.1%} complete)")
            
    result = pd.DataFrame(result)
    return result
    

        
if __name__ == "__main__":
    adata_path = sys.argv[1]
    annotations_path = sys.argv[2]
    output_path = sys.argv[3]
    
    # Load the data
    print(f"\nLoading data from {adata_path}...")
    adata = sc.read_h5ad(adata_path)
    
    # Load cell types
    print(f"\nLoading annotations from {annotations_path}...")
    df = ut.load_annotations(annotations_path)
    
    print(f"\nMerging annotations...")
    adata.obs = pd.merge(
        adata.obs,
        df,
        how='left',
        left_index=True,
        right_index=True,
    )
    
    print(f"\nTokenizing data...")
    result = tokenize(adata)

    # Save the data
    print(f"\nSaving data to {output_path}...")

    dataset = Dataset.from_pandas(result)
    dataset.save_to_disk(output_path)
    print("Done!")
    
    

 
    
    
    
    