#!/bin/bash
#SBATCH --job-name=ae_sgts_new
#SBATCH --output=/scratch/nenarokova/code/slurm_out/ae_sgts_new_%A_%a.out
#SBATCH --partition=high
#SBATCH --time=99-99:00:00
#SBATCH --array=1-85
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
##SBATCH --nodes=1
## --cpu_bind=v,threads

workdir="/scratch/nenarokova/euk/markers/ae/one_hit/01_07_23/no_idun/"
fasta_dir=$workdir"faa/"

linsi_dir=$workdir"linsi/"
trimmed_linsi_dir=$workdir"linsi_bmge/"

cd $fasta_dir

fasta=$(ls *.faa | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo $fasta

msa=$linsi_dir$fasta
trimmed_msa=$trimmed_linsi_dir$fasta

linsi --anysymbol $fasta > $msa
BMGE -i $msa -t "AA" -m BLOSUM30 -of $trimmed_msa

cd $trimmed_linsi_dir

iqtree2 -s $fasta -m LG+G -B 1000 -nt 1 -redo