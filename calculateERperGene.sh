#!/bin/bash

# Script to extract a gene from a multi fasta alignment file and calculate 
# the average number of substitutions per site per year for two timepoints (in that gene).

# sh calculateERperGene.sh CVR7347_withref.fasta sequence.gb S 20/4/2021 18/5/2021

alignment=$1
genbank=$2
gene=$3
date1=$4
date2=$5

python3 svartools.py extract -aln ${alignment} -gb ${genbank} -g ${gene} 

for i in ${gene}*.fasta; do
	echo "Opening file $i"
	python3 calculateER.py -f $i -d1 ${date1} -d2 ${date2}
done

rm ${gene}*.fasta