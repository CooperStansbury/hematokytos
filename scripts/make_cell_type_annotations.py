import sys
import os
import pandas as pd
import numpy as np
import yaml


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


def process_cell_type_files(file_list, cell_type_columns):
    """
    Processes a list of CSV files, extracts cell type information, and combines them into a single DataFrame.

    Args:
        file_list: A list of file paths to CSV files containing cell type data.
        cell_type_columns: A dictionary mapping original column names to standardized names.

    Returns:
        A DataFrame containing 'obs_index', 'cell_type', and 'dataset' columns.
    """
    result = []

    for fpath in file_list:
        basename = os.path.basename(fpath)
        sample_id = basename.split(".")[0]

        df = pd.read_csv(fpath)
        df = df.rename(columns=cell_type_columns)

        # If 'cell_type' doesn't exist, create it using sample_id
        if 'cell_type' not in df.columns:
            df['cell_type'] = sample_id

        df = df[['obs_index', 'cell_type']]
        df['dataset'] = sample_id

        result.append(df)

    return pd.concat(result)




        
if __name__ == "__main__":
    output_path = sys.argv[1]
    cell_column_path = sys.argv[2]
    cell_types_path = sys.argv[3]
    file_list = sys.argv[4:]

    # Load configuration files
    print("\nLoading cell column mapping...")
    cell_type_columns = load_yaml(cell_column_path)
    print("Cell column mapping loaded successfully!")

    print("\nLoading cell types...")
    cell_types = load_yaml(cell_types_path)
    print("Cell types loaded successfully!")

    # Process cell type files
    print("\nProcessing cell type files...")
    df = process_cell_type_files(file_list, cell_type_columns)
    print(f"Processed {len(file_list)} files.")

    # Standardize cell types
    print("\nStandardizing cell types...")
    df['standard_cell_type'] = df['cell_type'].map(cell_types)
    df['standard_cell_type'] = df['standard_cell_type'].fillna('NA')

    # Save the data
    print(f"\nSaving standardized cell types to {output_path}...")
    df.to_csv(output_path, index=False)
    print("Data saved successfully!")
    print("Done!")

 
    
    
    
    