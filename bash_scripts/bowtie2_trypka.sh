#!/bin/bash

bw2_dir='/home/nenarokova/tools/bowtie2-2.2.9/'
base_name='/media/4TB1/blastocrithidia/mapping/jac_bowtie2_RNA/jac_RNA_bw2'
ref="/media/4TB1/blastocrithidia/transcriptome_assembly/trinity_denovo/jac_default/Trinity.fasta"
$bw2_dir'bowtie2-build' --threads 32 $ref $base_name

p1_1="/media/4TB1/blastocrithidia/jac_rna_reads/trimmed/jac_trimmed_1.fq"
p1_2="/media/4TB1/blastocrithidia/jac_rna_reads/trimmed/jac_trimmed_1.fq"

alignment=$base_name".sam"
report=$base_name".txt"
unmapped_unpaired=$base_name"_unmapped_unpaired.fq"
unmapped_paired=$base_name"_unmapped_paired.fq"

$bw2_dir'bowtie2' --very-sensitive -p 30 -x $base_name -1 $p1_1 -2 $p1_2 --un-gz $unmapped_unpaired --un-conc-gz $unmapped_paired -S $alignment 2> $report

samfile=$alignment
bamfile=$base_name"_unsorted.bam"
sorted=$base_name"_sorted"
sorted_file=$sorted".bam"

samtools view -bS $samfile > $bamfile -@ 20
samtools sort -o $sorted_file -@ 20 $bamfile
samtools index -b $sorted_file