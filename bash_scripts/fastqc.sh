#!/bin/sh

read_dir='/home/kika/diplonema/reads/merged_qtrimmed_vstrict/'
out_dir='/home/kika/diplonema/reads/merged_qtrimmed_vstrict/fastqc/'

/home/kika/tools/FastQC/fastqc -o $out_dir $read_dir'YPF1604_adapter_trimmed_merged.fq'
/home/kika/tools/FastQC/fastqc -o $out_dir $read_dir'YPF1604_adapter_trimmed_unmerged_1.fq'
/home/kika/tools/FastQC/fastqc -o $out_dir $read_dir'YPF1604_adapter_trimmed_unmerged_2.fq'