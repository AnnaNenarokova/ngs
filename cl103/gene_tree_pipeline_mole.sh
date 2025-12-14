#!/bin/bash
#SBATCH --job-name=bac_markers_trees
#SBATCH --output=/home/users/nenarokova/slurm_output/bac_markers_trees_%A_%a.out
#SBATCH --time=99-99:00:00
#SBATCH --array=1-114
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB

set -euo pipefail

module purge
module load mafft
module load BMGE-1.12
module load iqtree-2.1.0

workdir="/scratch/nenarokova/bacteria/markers/"
fasta_dir="${workdir}/fasta/"
linsi_dir="${workdir}/linsi/"
trimmed_dir="${workdir}/linsi_bmge/"
tree_dir="${workdir}/treefiles_LGG/"
final_tree_dir="${workdir}/trees_LGG/"

mkdir -p "$linsi_dir" "$trimmed_dir" "$tree_dir" "$final_tree_dir"

echo "[$(date)] Job started on node $(hostname)"
echo "[$(date)] Selecting FASTA for array task ${SLURM_ARRAY_TASK_ID}"

cd "$fasta_dir"

fasta=$(ls *.fasta *.faa 2>/dev/null | sed -n "${SLURM_ARRAY_TASK_ID}p" || true)

if [[ -z "${fasta:-}" ]]; then echo "[$(date)] ERROR: No FASTA found for task ${SLURM_ARRAY_TASK_ID}"; exit 1; fi

echo "[$(date)] FASTA selected: $fasta"

msa="${linsi_dir}/${fasta}.aln"
trimmed_msa="${trimmed_dir}/${fasta}.bmge.aln"

echo "[$(date)] START MAFFT L-INS-i alignment"
mafft --maxiterate 1000 --localpair --anysymbol "$fasta" > "$msa"
echo "[$(date)] END MAFFT alignment: $msa"

echo "[$(date)] START BMGE trimming"
java -jar /home/software/BMGE-1.12/BMGE.jar -i "$msa" -t AA -m BLOSUM30 -of "$trimmed_msa"
echo "[$(date)] END BMGE trimming: $trimmed_msa"

echo "[$(date)] Copying trimmed alignment to IQ-TREE directory"
cp "$trimmed_msa" "$tree_dir/"

cd "$tree_dir"

echo "[$(date)] START IQ-TREE inference"
iqtree2 -s "$(basename "$trimmed_msa")" -m LG+G -B 1000 -nt 1 -redo
echo "[$(date)] END IQ-TREE inference"

echo "[$(date)] Copying final tree to final tree directory"
cp "$(basename "$trimmed_msa").treefile" "$final_tree_dir/"

echo "[$(date)] Finished marker: $fasta"
echo "[$(date)] Job completed"