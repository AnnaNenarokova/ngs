#!/bin/bash
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=60

folder='/home/nenarokova/contaminants/trimmed_reads'
cd $folder
alignment='/home/nenarokova/contaminants/genomes/test.sam'

bt2_base='/home/nenarokova/contaminants/genomes/seymouri'

/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 H10-3_ACAGTG_L007_001_ad_q20_l50_paired_out_fw.fastq -2 H10-3_ACAGTG_L007_001_ad_q20_l50_paired_out_rv.fastq -U H10-3_ACAGTG_L007_001_ad_q20_l50_unpaired_out_fw.fastq,H10-3_ACAGTG_L007_001_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/h10_seymouri_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_seymouri_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 Leptomonas_seymouri_ad_q20_l50_paired_out_fw.fastq -2 Leptomonas_seymouri_ad_q20_l50_paired_out_rv.fastq -U Leptomonas_seymouri_ad_q20_l50_unpaired_out_fw.fastq,Leptomonas_seymouri_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/seymouri_seymouri_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 B08_TAGCTT_L001_ad_q20_l50_paired_out_fw.fastq -2 B08_TAGCTT_L001_ad_q20_l50_paired_out_rv.fastq -U B08_TAGCTT_L001_ad_q20_l50_unpaired_out_fw.fastq,B08_TAGCTT_L001_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/b08_seymouri_vvfast.txt

bt2_base='/home/nenarokova/contaminants/genomes/leptomonas'

/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_h10_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 Leptomonas_seymouri_ad_q20_l50_paired_out_fw.fastq -2 Leptomonas_seymouri_ad_q20_l50_paired_out_rv.fastq -U Leptomonas_seymouri_ad_q20_l50_unpaired_out_fw.fastq,Leptomonas_seymouri_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/seymouri_h10_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 B08_TAGCTT_L001_ad_q20_l50_paired_out_fw.fastq -2 B08_TAGCTT_L001_ad_q20_l50_paired_out_rv.fastq -U B08_TAGCTT_L001_ad_q20_l50_unpaired_out_fw.fastq,B08_TAGCTT_L001_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/b08_h10_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 H10-3_ACAGTG_L007_001_ad_q20_l50_paired_out_fw.fastq -2 H10-3_ACAGTG_L007_001_ad_q20_l50_paired_out_rv.fastq -U H10-3_ACAGTG_L007_001_ad_q20_l50_unpaired_out_fw.fastq,H10-3_ACAGTG_L007_001_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/h10_h10_vvfast.txt

bt2_base='/home/nenarokova/contaminants/genomes/blechomonas'

/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 c_thermophila_dna_ad_q20_l50_paired_out_fw.fastq -2 c_thermophila_dna_ad_q20_l50_paired_out_rv.fastq -U c_thermophila_dna_ad_q20_l50_unpaired_out_fw.fastq,c_thermophila_dna_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/thermophila_b08_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 Leptomonas_seymouri_ad_q20_l50_paired_out_fw.fastq -2 Leptomonas_seymouri_ad_q20_l50_paired_out_rv.fastq -U Leptomonas_seymouri_ad_q20_l50_unpaired_out_fw.fastq,Leptomonas_seymouri_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/seymouri_b08_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 B08_TAGCTT_L001_ad_q20_l50_paired_out_fw.fastq -2 B08_TAGCTT_L001_ad_q20_l50_paired_out_rv.fastq -U B08_TAGCTT_L001_ad_q20_l50_unpaired_out_fw.fastq,B08_TAGCTT_L001_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/b08_b08_vvfast.txt
/home/nenarokova/tools/bowtie2-2.2.9/bowtie2 -D 5 -R 1 -N 0 -L 32 -p 60 -x $bt2_base -1 H10-3_ACAGTG_L007_001_ad_q20_l50_paired_out_fw.fastq -2 H10-3_ACAGTG_L007_001_ad_q20_l50_paired_out_rv.fastq -U H10-3_ACAGTG_L007_001_ad_q20_l50_unpaired_out_fw.fastq,H10-3_ACAGTG_L007_001_ad_q20_l50_unpaired_out_rv.fastq -S $alignment 2> /home/nenarokova/contaminants/bw2_stats/h10_b08_vvfast.txt
