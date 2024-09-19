rule get_gtf:
    input:
        config['gtf_path']
    output:
        OUTPUT + "reference/annotations.gtf"
    shell:
        """cp {input} {output}"""
        

rule make_gene_table:
    input:
        OUTPUT + "reference/annotations.gtf"
    output:
        OUTPUT + "reference/gene_table.csv"
    conda:
        "bioinf"
    shell:
        """python scripts/make_gene_table.py {input} {output}"""
        
        
rule make_token_map:
    input:
        gene_table=OUTPUT + "reference/gene_table.csv",
        tokens=OUTPUT + "reference/token_dictionary.pkl",
    output:
        OUTPUT + "reference/token_map.csv"
    conda:
        "bioinf"
    shell:
        """python scripts/make_token_map.py {input.gene_table} {input.tokens} {output}"""
        

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
        train=OUTPUT + "pretrained_model/training_args.bin",
        config=OUTPUT + "pretrained_model/config.json",
        model=OUTPUT + "pretrained_model/pytorch_model.bin",
    params:
        output_dir=OUTPUT + "pretrained_model/"
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
        