from datetime import datetime
import pandas as pd
import yaml
from pathlib import Path
import re
import os
import sys
from tabulate import tabulate


BASE_DIR = Path(workflow.basedir)
configfile: str(BASE_DIR) + "/config/config.yaml"

# big picture variables
OUTPUT = config['output_path']

# load in images paths
input_path = os.path.abspath(config['inputs'])
input_df = pd.read_csv(input_path, comment="#")
samples = input_df['sample_id'].to_list()

# get new path names
output_paths = [OUTPUT + "raw_anndata/" + x + '.h5ad' for x in samples]

# print statements
print("\n----- CONFIG VALUES -----")
for key, value in config.items():
    print(f"{key}: {value}")
    
    
print("\n----- INPUT VALUES -----")
print(
    tabulate(
        input_df, 
        headers='keys', 
        tablefmt='psql',
        showindex=False,
    )
)


rule all:
    input:
        expand(OUTPUT + "raw_anndata/{sid}.h5ad", sid=samples),
        OUTPUT + "model/training_args.bin",
        OUTPUT + "reference/token_dictionary.pkl",
        OUTPUT + "reference/gene_median_dictionary.pkl",
        
   
rule get_raw_anndata:
    input:
        input_df['file_path'].to_list()
    output:
        output_paths
    run:
        from shutil import copyfile
        for i, refPath in enumerate(input):

            outPath = output[i]
            copyfile(refPath, outPath)
            
            
rule get_model_weights:
    input:
        config['model_path']
    output:
        train=OUTPUT + "model/training_args.bin",
        config=OUTPUT + "model/config.json",
        model=OUTPUT + "model/pytorch_model.bin",
    params:
        output_dir=OUTPUT + "model/"
    shell:
        """cp -r {input}/* {params.output_dir}"""
            
        
        
rule get_token_dict:
    input:
        config['token_dict_path']
    output:
        OUTPUT + "reference/token_dictionary.pkl"
    shell:
        """cp {input} {output}"""
        
        
         
rule get_gene_lengths:
    input:
        config['gene_norm_path']
    output:
        OUTPUT + "reference/gene_median_dictionary.pkl"
    shell:
        """cp {input} {output}"""