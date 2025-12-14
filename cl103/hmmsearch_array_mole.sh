#!/bin/bash
#SBATCH --job-name=gtdb_markers_hmm
#SBATCH --output=/home/users/nenarokova/slurm_output/gtdb_markers_hmm_%A_%a.out
#SBATCH --time=7-12:00:00
#SBATCH --array=1-133
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB

module load hmmer3.3

proteome_dir="/home/users/nenarokova/daria/cl103/phylogenetics/proteomes_bacteria/"
hmm_db="/home/users/nenarokova/daria/gtdb/release226/markers/tigrfam/tigrfam.hmm"
out_dir="/home/users/nenarokova/daria/cl103/phylogenetics/hmms/hmm_markers_results/"

mkdir -p "$out_dir"

cd "$proteome_dir"
proteome=$(ls *.fasta | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$proteome" ]]; then
    echo "No proteome found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

proteome_path="${proteome_dir}${proteome}"
tbl_out="${out_dir}${proteome}.tbl"
domtbl_out="${out_dir}${proteome}.domtbl"
log_out="${out_dir}${proteome}.hmmsearch.log"

echo "Proteome:  $proteome_path"
echo "HMM DB:    $hmm_db"
echo "tblout:    $tbl_out"
echo "domtblout: $domtbl_out"

hmmsearch --noali --cut_ga --tblout "$tbl_out" --domtblout "$domtbl_out" "$hmm_db" "$proteome_path" > "$log_out"