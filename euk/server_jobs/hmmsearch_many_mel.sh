#!/bin/bash
#SBATCH --job-name=hmmsearch
#SBATCH --output=/scratch/nenarokova/code/slurm_out/hmmsearch_%A_%a.out
#SBATCH --partition=high
#SBATCH --time=7-12:00:00
#SBATCH --array=1-85
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
##SBATCH --nodes=1
## --cpu_bind=v,threads

hmmsearch="/mnt/alvarium2pool/scratch/nenarokova/tools/hmmer-3.3.2/bin/hmmsearch"

hmm_dir="/scratch/nenarokova/euk/markers/common_results/mono_euk_sets/set3/hmm/"

subject_dir=""
hmm_results_dir="/scratch/nenarokova/euk/markers/ae/one_hit/hmmresults_monoeuks_10_06/"

cd $fasta_dir
fasta=$(ls *.faa | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo $fasta

msa=$msa_dir$fasta
hmm_file=$hmm_dir$fasta".hmm"

e_threshold="0.0000001"

cd $subject_dir
for subject in *.fasta
do
	subject_path=$subject_dir$subject
	result=$hmm_results_dir$subject$fasta".txt"
	echo $result
	hmmsearch -E $e_threshold --tblout $result $hmm_file $subject_path
done