#!/bin/bash
fasta="/Users/vl18625/work/novymonas_ncbi/nesm_ncbi.fasta"
gff="/Users/vl18625/work/novymonas_ncbi/nesm_edited_for_ncbi.gff"
template="/Users/vl18625/work/novymonas_ncbi/template.sbt"
out_sqn="/Users/vl18625/work/novymonas_ncbi/n_esmeraldas.sqn"
table2asn="/Users/vl18625/work/tools/mac.table2asn"
error_out="/Users/vl18625/work/novymonas_ncbi/t2asn_errors_10_01_22.txt"
$table2asn -M n -J -c w -euk -t $template -gaps-min 10 -l paired-ends -i $fasta -f $gff -o $out_sqn -locus-tag-prefix -augustus-fix -j "[organism=Novymonas esmeraldas][strain=E262AT.01]" -Z 2> $error_out


