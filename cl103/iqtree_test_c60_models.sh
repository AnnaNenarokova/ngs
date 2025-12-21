#!/bin/bash
#SBATCH --job-name=iqtree_C60_MF
#SBATCH --output=/home/users/nenarokova/slurm_output/iqtree_C60_MF_%A.out
#SBATCH --time=7-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G

set -euo pipefail

module load iqtree-2.1.0

workdir="/home/users/nenarokova/daria/cl103/phylogenetics/152_gtdb_markers/supermatrix/model_tester"
msa="/home/users/nenarokova/daria/cl103/phylogenetics/152_gtdb_markers/supermatrix/114_gtdb_markers_133_proteomes.fasta"

cd $workdir
iqtree2 -s "$msa" -T 8 -m MF -mset LG -madd C60 -mrate G --score-diff all -mem 250G -redo -pre C60_modeltest
