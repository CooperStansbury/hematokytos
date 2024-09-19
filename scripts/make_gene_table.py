import pandas as pd
import pyranges as pr
import pickle
import os
import sys
import time


def load_gtf(gtf_path):
    """
    Loads gene information from a GTF file into a Pandas DataFrame.

    Args:
        gtf_path: The path to the GTF file.

    Returns:
        A Pandas DataFrame containing the following columns for each gene:
            - gene_name
            - gene_id
            - gene_biotype
            - Chromosome
            - Start
            - End
    """

    keep_columns = ['gene_name', 'gene_id', 'gene_biotype', 'Chromosome', 'Start', 'End']
    gdf = pr.read_gtf(gtf_path).as_df()
    gdf = gdf[(gdf['Feature'] == 'gene') & gdf['gene_id'].notna()][keep_columns]
    return gdf.drop_duplicates().reset_index(drop=True)


        
if __name__ == "__main__":
    gtf_path = sys.argv[1]
    output_path = sys.argv[2]

    # Load the data
    print(f"\nLoading data from {gtf_path}...")
    start_time = time.time()
    gdf = load_gtf(gtf_path)
    end_time = time.time()
    print(f"Data loading completed in {end_time - start_time:.2f} seconds.")

    # Save the data
    print(f"\nSaving gene table to {output_path}...")
    start_time = time.time()
    gdf.to_csv(output_path, index=False)
    end_time = time.time()
    print(f"Saving completed in {end_time - start_time:.2f} seconds.")
    print("Done!")
 
    
    
    
    