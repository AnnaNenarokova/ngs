#!/bin/bash
#SBATCH --job-name=iqtree_C60_compare
#SBATCH --output=/home/users/nenarokova/slurm_output/iqtree_C60_compare_%A.out
#SBATCH --time=7-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G

set -euo pipefail

module load iqtree-2.1.0

msa="/home/users/nenarokova/daria/cl103/phylogenetics/152_gtdb_markers/supermatrix/114_gtdb_markers_133_proteomes.fasta"

iqtree2 -s "$msa" -T 8 -m LG+C60+G4,LG+C60+G8,LG+C60+F+G4,LG+C60+F+G8 --score-diff all -mem 250G -redo -pre C60_model_compare