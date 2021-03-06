#!/bin/bash
#PBS -l nodes=1:ppn=64

align_dir='/home/nenarokova/genomes/novymonas/pandoraea_phylogeny/all_OGs_1prot_short_no_Limno_alignment/'
raxml='/home/nenarokova/tools/raxml/raxml_sse/standard-RAxML/raxmlHPC-SSE3'

cd $align_dir
for f in *".fa.phy.trimal_auto1"
do
    $raxml -m PROTGAMMAGTR -p 12345 -s $f -n $f"2.tre" > $f"raxml_out.txt" &
done
