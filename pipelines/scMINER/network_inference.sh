#!/bin/bash

#SBATCH --job-name=scMINER
#SBATCH --account=indikar1
#SBATCH --partition=largemem,standard
#SBATCH --mail-user=cstansbu@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH --nodes=1                     
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=8

set -euo pipefail

base_dir="/nfs/turbo/umms-indikar/shared/projects/HSC/pipeline_outputs/scMINER/test/SJARACNe"
p_consensus="0.001"
p_bootstrap="0.0001"

echo "=============================================="
echo "Starting SJARACNe runs for C1–C5"
echo "Base directory: $base_dir"
echo "Consensus (p-value): $p_consensus"
echo "Bootstrap (p-value): $p_bootstrap"
echo "----------------------------------------------"

for cdir in "$base_dir"/C[1-5]; do
    cname=$(basename "$cdir")

    echo ""
    echo ">>> Processing cluster: $cname"

    exp_file=$(find "$cdir" -maxdepth 1 -name "${cname}.*.exp.txt")
    tf_file=$(find "$cdir/TF" -maxdepth 1 -name "${cname}.*.tf.txt")
    out_dir="$cdir/output"
    tmp_dir="/scratch/indikar_root/indikar1/cstansbu/tmp"

    echo "  - Expression file: $exp_file"
    echo "  - TF file        : $tf_file"
    echo "  - Output dir     : $out_dir"
    echo "  - Temp dir       : $tmp_dir"

    mkdir -p "$out_dir"

    echo "  > Running sjaracne..."
    sjaracne local \
      --exp-file "$exp_file" \
      --hub-genes "$tf_file" \
      --output-dir "$out_dir" \
      --tmpdir-prefix "$tmp_dir" \
      -pc "$p_consensus" \
      -pb "$p_bootstrap"

    echo "  ✓ Finished $cname"
done

echo ""
echo "=============================================="
echo "All SJARACNe runs complete."
