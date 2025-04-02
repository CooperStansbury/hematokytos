#!/bin/bash

#SBATCH --job-name=hsc_epi2me 
#SBATCH --account=indikar99
#SBATCH --partition=standard
#SBATCH --mail-user=cstansbu@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=150G
#SBATCH --time=36:00:00
#SBATCH --nodes=1                     
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=16 

conda activate hsc_epi2me
python pipeline_runner.py --force 