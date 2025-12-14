#!/bin/bash
#SBATCH --job-name=write_152_marker_fastas
#SBATCH --output=/home/users/nenarokova/slurm_output/write_152_marker_fastas_%j.out
#SBATCH --time=99:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Load Python
PYTHON="/usr/bin/python"
SCRIPT="/home/users/nenarokova/ngs/cl103/write_marker_fastas_hmmresults.py"

TBL_DIR="/home/users/nenarokova/daria/cl103/phylogenetics/hmms/hmm_markers_results/"
PROTEOMES_DIR="/home/users/nenarokova/daria/cl103/phylogenetics/proteomes_bacteria/"
OUT_DIR="/home/users/nenarokova/daria/cl103/phylogenetics/152_gtdb_markers/fasta/"

echo "Job started on $(hostname)"
echo "Start time: $(date)"

$PYTHON "$SCRIPT" \
  --tbl_dir "$TBL_DIR" \
  --proteomes_dir "$PROTEOMES_DIR" \
  --out_dir "$OUT_DIR" \
  --resolve_conflicts

echo "End time: $(date)"
echo "Job finished"
