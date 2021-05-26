# calculateER
Calculate the rate of evolution from two sequences between two time points.

#### Usage:
```
usage: calculateER.py [-h] [-f FASTA] [-d1 DATE1] [-d2 DATE2]

Calculate the rate of evolution between two time points.

optional arguments:

  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Path to aligned multifasta file with the two consensus sequences
  -d1 DATE1, --date1 DATE1
                        Respective date of the first sequence in the format of: dd/mm/yyyy
  -d2 DATE2, --date2 DATE2
                        Respective date of the second sequence in the format of: dd/mm/yyyy
```

#### Input:

* The -f input should be a multi-fasta file with the two aligned consensus sequences of interest. 

Example:

\>sequenceID-001 description \
AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT \
ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG \
\>sequenceID-002 description \
CAGTAAAGAGTGGATGTAAGAACCGTCCGATCTACCAGATGTGATAGAGGTTGCCAGTAC \
AAAAATTGCATAATAATTGATTAATCCTTTAATATTGTTTAGAATA

* The date input (-d1 and -d2) should have the following format:
dd/mm/yyyy

Example:
6/5/2021


# calculateER
This pipeline uses [svartools](https://github.com/iouliadt/svartools) to extract a gene specified by the user and then calculates he average number of substitutions per site per year using calculateER.

```
usage: calculateERperGene.sh multifasta_withref.fasta genbankFile.gb gene date1 date2

```

#### Input:

* The first entry of the multifasta alignment must be the reference sequence.

Example:

\>referenceID description \
AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT \
ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG \
\>sequenceID-001 description \
AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT \
ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG \
\>sequenceID-002 description \
CAGTAAAGAGTGGATGTAAGAACCGTCCGATCTACCAGATGTGATAGAGGTTGCCAGTAC \
AAAAATTGCATAATAATTGATTAATCCTTTAATATTGTTTAGAATA