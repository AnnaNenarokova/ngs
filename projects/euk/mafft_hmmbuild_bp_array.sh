#!/bin/bash

module load apps/mafft/7.429 apps/hmmer/3.3.2

bmge="/sw/apps/BMGE-1.12/BMGE.jar"

fasta_dir="/user/work/vl18625/euk/nina_markers/faa/"
msa_dir="/user/work/vl18625/euk/nina_markers/msa/"
hmm_dir="/user/work/vl18625/euk/nina_markers/hmm_profiles/"

cd $fasta_dir
fasta=$(ls *.faa | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo $fasta

msa=$msa_dir$fasta
hmm_file=$hmm_dir$fasta".hmm"

mafft --anysymbol $fasta > $msa
hmmbuild $hmm_file $msa