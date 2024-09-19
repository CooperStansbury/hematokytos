import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import torch
import json
from geneformer import Classifier
from datasets import Dataset, load_from_disk
from datasets import load_dataset
from geneformer import EmbExtractor

# local imports
import utils as ut

sc.settings.verbosity = 3  

def print_nested_params(params, indent=1):
    """Recursively prints nested parameters with indentation."""
    for key, value in params.items():
        if isinstance(value, dict):
            print("  " * indent + f"{key}:")
            print_nested_params(value, indent + 1)
        else:
            print("  " * indent + f"{key}: {value}")
            
            
def filter_labels(data, params):
    """Filters labels based on exclusions from params and prints the results.

    Args:
        data: The dataset containing labels.
        params: A dictionary containing 'label_column' and 'exclude_labels'.
    """

    label_column = params['label_column']
    labels = data.unique(label_column)
    exclusions = params['exclude_labels']
    filter_labels = [x for x in labels if x not in exclusions]

    filter_data_dict = {
        label_column: filter_labels, 
    }

    # Print filtered labels and exclusions 
    print("\n--- Filtering Labels ---")
    print(f"Unique Labels: {labels}")
    print(f"Excluded Labels: {exclusions}")
    print(f"Filtered Labels: {filter_labels}\n")

    return filter_data_dict  # Return the filtered data dictionary


        
if __name__ == "__main__":
    data_path = sys.argv[1]
    model_path = sys.argv[2]
    params_path = sys.argv[3]
    output_path = sys.argv[4]

    # some CUDA set up
    torch.cuda.empty_cache()
    CORES = os.cpu_count()
    GPUS = torch.cuda.device_count()

    # load config
    params = ut.load_yaml(params_path)

    # Nicely formatted print statements
    print("\n--- Configuration ---")
    print(f"Data Path: {data_path}")
    print(f"Model Path: {model_path}")
    print(f"Params Path: {params_path}")
    print(f"Output Path: {output_path}")
    print(f"Available Cores: {CORES}")
    print(f"Available GPUs: {GPUS}")
    print("--- Loaded Parameters ---")
    print("\n--- Loaded Parameters ---")
    print_nested_params(params)
    print("----------------------\n")

     # Load the data
    print(f"\n--- Loading data from {data_path} ---") 
    data = load_from_disk(data_path)
    print(f"Data loaded successfully. Shape: {data.shape}")  
    
    # filter the labels
    filter_data_dict = filter_labels(data, params)
    

    training_args = params['training_args']
    sample_size = None

    cc = Classifier(
        classifier = "cell",
        cell_state_dict = cell_state_dict,
        training_args = training_args,
        filter_data=filter_data_dict,
        max_ncells = None,
        freeze_layers = 2,
        num_crossval_splits = 1,
        forward_batch_size = 200,
        nproc = CORES,
        ngpu = GPUS,
    )
    
    
    # set up the 

#     label_column = params['label_column']
#     labels = data.unique(label_column)
#     exclusions = params['exclude_labels']
#     filter_labels = [x for x in labels if x not in exclusions]

#     filter_data_dict = {
#         label_column: filter_labels, 
#     }

#     # Print filtered labels and exclusions 
#     print("\n--- Filtering Labels ---")
#     print(f"Unique Labels: {labels}")
#     print(f"Excluded Labels: {exclusions}")
#     print(f"Filtered Labels: {filter_labels}\n")
    
    
    
    
#     # Extract the observation metadata
#     df = adata.obs.copy()
#     df = df.reset_index(drop=False, names='obs_index')

#     # Save the data
#     print(f"\nSaving metadata to {output_path}...")
#     df.to_csv(output_path, index=False)
#     print("Done!")
    
    

 
    
    
    
    