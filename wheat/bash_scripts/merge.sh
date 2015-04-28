#!/bin/bash
cd /mnt/results/nenarokova/wheat/R/sum_fastq_re/sorted
mkdir merged_alignments
for l in A B C D E F G H
do
	out_bam='./merged_alignments/'$l$PBS_ARRAYID'.bam'
	in_bam1='./'$l$PBS_ARRAYID'/new_assembly_nbs_lrr_ids.bam'
	in_bam2='./'$l'0'$PBS_ARRAYID'/new_assembly_nbs_lrr_ids.bam'
	echo $out_bam
	echo $in_bam1
	echo $in_bam2
	# samtools merge $out_bam $in_bam1 $in_bam
done