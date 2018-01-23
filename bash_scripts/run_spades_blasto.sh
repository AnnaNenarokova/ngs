#!/bin/bash

read_dir='/home/kika/diplonema/reads/trimmed/'
pe1_1=$read_dir'YPF1604_trimmed_1.fq.gz'
pe1_2=$read_dir'YPF1604_trimmed_1.fq.gz'

outdir='/home/kika/diplonema/genome_assembly/1604/'
report=$outdir'spades_report.txt'

/home/kika/tools/SPAdes-3.11.1-Linux/bin/spades.py --pe1-1 $pe1_1 --pe1-2 $pe1_2 --careful -t 30 -o $outdir 2> $report
