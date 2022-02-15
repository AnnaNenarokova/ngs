#!/bin/bash

threads=30

fasta="/media/4TB1/blastocrithidia/genome_assembly/bnonstop_corrected_assembly.fasta"
bt2_base="/media/4TB1/blastocrithidia/mapping/bw2_indexes/bnonstop_corrected_bw2"

file_path="/media/4TB1/blastocrithidia/mapping/p57_ill_bw2_tRNA_DNA/new_tRNA_ill_p57"
reads="/media/4TB1/blastocrithidia/reads/tRNAs/trimmed/P57-cyto_trimmed.fq.gz"

aligned=$file_path"_aligned.fq.gz"

alignment=$file_path".sam"
report=$file_path".txt"

bowtie2-build --threads $threads $fasta $bt2_base
bowtie2 --very-sensitive -p $threads -x $bt2_base -U $reads --al $aligned -S $alignment 2> $report

samfile=$alignment
unsorted_bam=$file_path"_unsorted.bam"
sorted_bam=$file_path"_sorted.bam"

samtools view -bS -@ $threads $samfile > $unsorted_bam
samtools sort -o $sorted_bam -@ $threads $unsorted_bam
samtools index -b $sorted_bam

