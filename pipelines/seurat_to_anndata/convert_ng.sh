#!/bin/bash

# Set input and output paths
inpath="/nfs/turbo/umms-indikar/shared/projects/HSC/data/datasets/ng_2024/iHSC.rds"
outpath="/nfs/turbo/umms-indikar/shared/projects/HSC/data/datasets/ng_2024/iHSC.h5ad"

# Run the R conversion script with arguments
Rscript convert.r "$inpath" "$outpath"