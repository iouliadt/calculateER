#!/bin/bash

# Script to extract gene from multi fasta alignment file and
# find the average number of substitutions per site per year for two timepoints.

# sh calculateERperGene.sh CVR7347_withref.fasta sequence.gb S 20/4/2021 18/5/2021

alignment=$1
genbank=$2
gene=$3
date1=$4
date2=$5

python3 svartools.py extract -aln ${alignment} -gb ${genbank} -g ${gene} 

for i in ${gene}*.fasta; do
	echo "Processing $i"
	python3 calculateER.py -f $i -d1 ${date1} -d2 ${date2}
done

rm ${gene}*.fasta