import pandas as pd
import pyranges as pr
import pickle
import os
import sys
import time


def load_gene_tokenization(token_dictionary_file):
    """
    Loads gene tokenization data from a pickle file into a DataFrame.

    Args:
        token_dictionary_file (str): Path to the pickle file containing the gene-token dictionary.

    Returns:
        pd.DataFrame: A DataFrame with 'Ensembl ID' and 'Token' columns.
    """

    with open(token_dictionary_file, "rb") as f:
        gene_token_dict = pickle.load(f)

    return pd.DataFrame(list(gene_token_dict.items()), columns=['gene_id', 'token_id'])

        
if __name__ == "__main__":
    gene_table_path = sys.argv[1]
    token_path = sys.argv[2]
    output_path = sys.argv[3]

    # Load the gene table
    print(f"\nLoading gene table from {gene_table_path}...")
    gdf = pd.read_csv(gene_table_path)
    print("Gene table loaded successfully!")

    # Load the token ids 
    print(f"\nLoading token data from {token_path}...")
    tokens = load_gene_tokenization(token_path)
    print("Token data loaded successfully!")

    # Merge
    print("\nMerging gene table and token data...")
    df = pd.merge(tokens, gdf, how='left')
    print("Merge completed!")

    # Save the data
    print(f"\nSaving token map to {output_path}...")
    df.to_csv(output_path, index=False)
    print("Token map saved successfully!")
    print("Done!")
 
    
    
    
    