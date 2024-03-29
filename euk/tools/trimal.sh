#!/bin/bash
trimal="/Users/vl18625/work/tools/trimAl/source/trimal"

al_dir="/Users/vl18625/work/euk/markers_euks/nina_markers/merged_arCOGs/linsi_bmge/"
cleaned_dir="/Users/vl18625/work/euk/markers_euks/nina_markers/merged_arCOGs/linsi_bmge_trimal_85/"

mkdir $cleaned_dir
cd $al_dir

for msa in *.faa
do
    cleaned_msa=$cleaned_dir$msa
    echo $msa
    $trimal -in $msa -out $cleaned_msa -resoverlap 0.85 -seqoverlap 85
done
