import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import torch
import json
import shutil
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
            
            
def manage_directory(directory_path):
    """
    Checks if a directory exists. If it does, removes it. Otherwise, creates it.

    Args:
        directory_path (str): The path to the directory to manage.
    """
    if os.path.exists(directory_path):
        if os.path.isdir(directory_path):
            shutil.rmtree(directory_path)  # Remove the directory and its contents
        else:
            os.remove(directory_path)  # Remove the file (if it's not a directory)

    os.makedirs(directory_path)  # Create the directory (including any necessary parent directories)
        
            
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
    print(f"Data Output Path: {output_path}")
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
    
    print(f"\n--- Setting up the Classifier ---") 
    
    cc = Classifier(
        classifier = "cell",
        cell_state_dict = params['cell_state_dict'],
        training_args = params['training_args'],
        filter_data=filter_data_dict,
        nproc = CORES,
        ngpu = GPUS,
        **params['classifier_args'],
    )
    
    # structuring output directory
    print(f"\n--- Building output directory: {output_path} ---") 
    manage_directory(output_path)
    print(os.listdir(output_path))
    
    print(f"\n--- Preparing the data ---") 
    torch.cuda.empty_cache()

    cc.prepare_data(
        input_data_file=data_path,
        output_directory=output_path,
        **params['data_preparation'],
    )
    
    print(f"\n--- Trainning GeneFormer ---") 
    torch.cuda.empty_cache()

    prefix = params['data_preparation']['output_prefix']
    root = f"{output_path}/{prefix}"

    all_metrics = cc.validate(
        model_directory=model_path,
        prepared_input_data_file=f"{root}_labeled_train.dataset",
        id_class_dict_file=f"{root}_id_class_dict.pkl",
        output_directory=output_path,
        save_eval_output=True,
        predict_trainer=True,
        n_hyperopt_trials=0,
        output_prefix=prefix,
    )
    
    print_nested_params(all_metrics)
    
    print("------------------------------------------")
    print("Done!")
    
    

 
    
    
    
    