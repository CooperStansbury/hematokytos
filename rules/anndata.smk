rule extract_obs:
    input:
        OUTPUT + "raw_anndata/{sid}.h5ad",
    output:
        OUTPUT + "obs/{sid}.csv"
    conda:
        "scanpy"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """python scripts/extract_obs.py {input} {output}"""
        

rule make_cell_type_annotations:
    input:
        obs=expand(OUTPUT + "obs/{sid}.csv", sid=samples),
        cell_columns=config['cell_column_map'],
        cell_types=config['cell_type_map'],
    output:
        OUTPUT + "annotation/cell_types.csv"
    conda:
        "bioinf"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """python scripts/make_cell_type_annotations.py {output} \
        {input.cell_columns} {input.cell_types} {input.obs}"""
    
        
rule process_anndata:
    input:
        anndata=OUTPUT + "raw_anndata/{sid}.h5ad",
        tokens=OUTPUT + "reference/token_map.csv",
    output:
        OUTPUT + "processed_anndata/{sid}.h5ad"
    conda:
        "scanpy"
    params:
        min_cells=config['min_cells'],
        min_genes=config['min_genes'],
        target_sum=config['target_sum'],
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    shell:
        """python scripts/process_anndata.py {input.anndata} {params.min_cells} \
        {params.min_genes} {params.target_sum} {input.tokens} {output}"""
        
        
rule get_qc_metrics:
    input:
        OUTPUT + "processed_anndata/{sid}.h5ad",
    output:
        obs=OUTPUT + "qc_metrics/{sid}.obs.csv",
        var=OUTPUT + "qc_metrics/{sid}.var.csv",
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    conda:
        "scanpy"
    shell:
        """python scripts/get_qc_metrics.py {input} {output.obs} {output.var}"""
        
        
rule merge_anndata:
    input:
        anndata=expand(OUTPUT + "processed_anndata/{sid}.h5ad", sid=samples),
        tokens=OUTPUT + "reference/token_map.csv",
    output:
        OUTPUT + "merged_anndata/merged_adata.h5ad"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    conda:
        "scanpy"
    shell:
        """python scripts/merge_anndata.py {output} {input.tokens} {input.anndata}"""    
        
       
rule combat_anndata:
    input:
        OUTPUT + "merged_anndata/merged_adata.h5ad"
    output:
        OUTPUT + "merged_anndata/combat_adata.h5ad"
    wildcard_constraints:
        sid='|'.join([re.escape(x) for x in set(samples)]),
    conda:
        "scanpy"
    shell:
        """python scripts/combat_anndata.py {input} {output}"""
        
        
rule tokenize_anndata:
    input:
        anndata=OUTPUT + "merged_anndata/{data}.h5ad",
        annotations=OUTPUT + "annotation/cell_types.csv"
    output:
         directory(OUTPUT + "tokenized_data/{data}.dataset")
    conda:
        "geneformer"
    shell:
        """python scripts/tokenize_anndata.py {input.anndata} \
        {input.annotations} {output}"""
               
        
        
        
        