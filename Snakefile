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


############### INPUT RULES ####################
include: "rules/gather.smk"
include: "rules/anndata.smk"


rule all:
    input:
        expand(OUTPUT + "raw_anndata/{sid}.h5ad", sid=samples),
        expand(OUTPUT + "processed_anndata/{sid}.h5ad", sid=samples),
        expand(OUTPUT + "obs/{sid}.csv", sid=samples),
        expand(OUTPUT + "qc_metrics/{sid}.obs.csv", sid=samples),
        OUTPUT + "pretrained_model/training_args.bin",
        OUTPUT + "reference/token_dictionary.pkl",
        OUTPUT + "reference/gene_median_dictionary.pkl",
        OUTPUT + "reference/annotations.gtf",
        OUTPUT + "reference/gene_table.csv",
        OUTPUT + "reference/token_map.csv",
        OUTPUT + "annotation/cell_types.csv",
        OUTPUT + "merged_anndata/merged_adata.h5ad",
        OUTPUT + "merged_anndata/combat_adata.h5ad",
        expand(OUTPUT + "tokenized_data/{data}.dataset", data=['merged_adata', 'combat_adata']),
        
        
rule gpu:
    input:
        expand(OUTPUT + "pretrained_embeddings/{data}.h5ad", data=['merged_adata', 'combat_adata']),
        
        
        
rule extract_pretrained_embeddings:
    input:
        data=OUTPUT + "tokenized_data/{data}.dataset",
        model=OUTPUT + "pretrained_model/"
    output:
        OUTPUT + "pretrained_embeddings/{data}.h5ad"
    conda:
        "geneformer"
    shell:
        """python scripts/extract_embeddings.py {input.data} \
        {input.model} {output} """
        