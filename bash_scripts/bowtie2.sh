#!/bin/bash
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=60

folder='/home/nenarokova/contaminants/trimmed_reads'
cd $folder
alignment='/home/nenarokova/contaminants/genomes/test.sam'

bt2_base='/home/nenarokova/contaminants/genomes/seymouri'

/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --very-fast -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_seymouri_vfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --fast -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_seymouri_fast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --sensitive -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_seymouri_sens.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --very-sensitive -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_seymouri_vsens.txt

bt2_base='/home/nenarokova/contaminants/genomes/leptomonas'

/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --very-fast -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_h10_vfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --fast -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_h10_fast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --sensitive -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_h10_sens.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --very-sensitive -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_h10_vsens.txt

bt2_base='/home/nenarokova/contaminants/genomes/blechomonas'

/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --very-fast -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_b08_vfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --fast -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_b08_fast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --sensitive -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_b08_sens.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 --very-sensitive -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_b08_vsens.txt
