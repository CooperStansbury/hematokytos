# scMINER Pipeline

Helper script to run gene regulatory network inference using the SJARACNe
implementation from scMINER.

## Contents
- `network_inference.sh` – SLURM batch script that loops over clusters and launches `sjaracne`.

## Usage
Edit the paths at the top of `network_inference.sh` to match your data
locations and submit it with `sbatch`:

```bash
sbatch network_inference.sh
```
