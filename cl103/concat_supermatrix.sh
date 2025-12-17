#!/bin/bash
#SBATCH --job-name=concat_supermatrix
#SBATCH --output=/scratch/nenarokova/code/slurm_out/concat_supermatrix_%j.out
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

module purge

msa_dir="/scratch/nenarokova/bacteria/markers/linsi_bmge/"
script="/scratch/nenarokova/bacteria/markers/scripts/concat_supermatrix.py"
out_fasta="/scratch/nenarokova/bacteria/markers/supermatrix/supermatrix.fasta"
out_part="/scratch/nenarokova/bacteria/markers/supermatrix/supermatrix.partitions"
msa_ext=".fasta.bmge.aln"
taxon_mode="before_pipe"

echo "[$(date)] Job started on $(hostname)"
echo "[$(date)] MSA dir: $msa_dir"
echo "[$(date)] Script: $script"
echo "[$(date)] Output FASTA: $out_fasta"
echo "[$(date)] Output partitions: $out_part"

python3 "$script" --msa_dir "$msa_dir" --msa_ext "$msa_ext" --taxon_mode "$taxon_mode" --out_fasta "$out_fasta" --partitions "$out_part"

echo "[$(date)] Job finished"