#!/bin/bash
#SBATCH --job-name=concat_supermatrix
#SBATCH --output=/home/users/nenarokova/slurm_output/concat_supermatrix_%j.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

module purge

script="/home/users/nenarokova/ngs/cl103/concat_supermatrix.py"
msa_dir="//home/users/nenarokova/daria/cl103/phylogenetics/152_gtdb_markers/linsi_bmge/"
out_fasta="/home/users/nenarokova/daria/cl103/phylogenetics/152_gtdb_markers/supermatrix/114_gtdb_markers_133_proteomes.fasta"
out_part="/home/users/nenarokova/daria/cl103/phylogenetics/152_gtdb_markers/supermatrix/114_gtdb_markers_133_proteomes.partitions"
msa_ext=".aln"
taxon_mode="before_pipe"

echo "[$(date)] Job started on $(hostname)"
echo "[$(date)] MSA dir: $msa_dir"
echo "[$(date)] Script: $script"
echo "[$(date)] Output FASTA: $out_fasta"
echo "[$(date)] Output partitions: $out_part"

python3 "$script" --msa_dir "$msa_dir" --msa_ext "$msa_ext" --taxon_mode "$taxon_mode" --out_fasta "$out_fasta" --partitions "$out_part"

echo "[$(date)] Job finished"