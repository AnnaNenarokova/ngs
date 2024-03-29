#!/bin/bash
#SBATCH --job-name=ae_81m_no_meta
#SBATCH --output=/scratch/nenarokova/code/slurm_out/ae_81m_no_meta_%A_%a.out
#SBATCH --partition=high
#SBATCH --time=99-99:00:00
#SBATCH --array=1-81
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
##SBATCH --nodes=1
## --cpu_bind=v,threads

workdir="/scratch/nenarokova/euk/markers/ae/one_hit/ae_no_meta_81_markers_11_03_24/"
fasta_dir=$workdir"faa/"

linsi_dir=$workdir"linsi/"
trimmed_linsi_dir=$workdir"linsi_bmge/"
lgg_dir=$workdir"treefiles_lgg/"

cd $fasta_dir

fasta=$(ls *.faa | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo $fasta

msa=$linsi_dir$fasta
trimmed_msa=$trimmed_linsi_dir$fasta
copy_trimmed_msa=$lgg_dir$fasta

linsi --anysymbol $fasta > $msa
BMGE -i $msa -t "AA" -m BLOSUM30 -of $trimmed_msa

cp $trimmed_msa $copy_trimmed_msa
cd $lgg_dir

iqtree2 -s $fasta -m LG+G -B 1000 -nt 1 -redo
