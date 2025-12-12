#!/bin/bash
#SBATCH --job-name=bact_markers_hmm
#SBATCH --output=/home/users/nenarokova/slurm_output/bact_marker_hmm_%A_%a.out
#SBATCH --time=7-12:00:00
#SBATCH --array=1-152
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
##SBATCH --nodes=1
## --cpu_bind=v,threads

module load hmmer3.3

subject_dir="/home/users/nenarokova/daria/cl103/phylogenetics/proteomes_bacteria/"
hmm_dir="/home/users/nenarokova/daria/cl103/phylogenetics/individual_hmms_tigr_gtdb/"
hmm_results_dir="/home/users/nenarokova/daria/cl103/phylogenetics/hmm_markers_results/"

cd $hmm_dir
hmm_file=$(ls *.hmm | sed -n ${SLURM_ARRAY_TASK_ID}p)
hmm_path=$hmm_dir$hmm_file
echo $hmm_path

e_threshold="0.0000001"

cd $subject_dir
for subject in *.fasta
do
	subject_path=$subject_dir$subject
	result=$hmm_results_dir$subject$hmm_file".txt"
	echo $result
	hmmsearch -E $e_threshold --tblout $result $hmm_path $subject_path
done