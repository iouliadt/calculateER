# calculateER
Calculate the rate of evolution from two sequences between two time points.

#### Usage:
```
usage: calculateER.py [-h] [-f FASTA] [-d1 DATE1] [-d2 DATE2]

Calculate the rate of evolution between two time points.

optional arguments:

  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Path to multifasta file with the two consensus sequences
  -d1 DATE1, --date1 DATE1
                        Respective date of the first sequence in the format of: dd/mm/yyyy
  -d2 DATE2, --date2 DATE2
                        Respective date of the second sequence in the format of: dd/mm/yyyy
```

#### Input:

* The -f input should be a multi-fasta file with the two consensus sequences of interest. 

Example:

\>sequenceID-001 description
AAGTAGGAATAATATCTTATCATTATAGATAAAAACCTTCTGAATTTGCTTAGTGTGTAT
ACGACTAGACATATATCAGCTCGCCGATTATTTGGATTATTCCCTG
\>sequenceID-002 description
CAGTAAAGAGTGGATGTAAGAACCGTCCGATCTACCAGATGTGATAGAGGTTGCCAGTAC
AAAAATTGCATAATAATTGATTAATCCTTTAATATTGTTTAGAATATATCCGTCAGATAA
TCCTAAAAATAACGATATGATGGCGGAAATCGTC

* The date input should have the following format:
dd/mm/yyyy

Example:
6/5/2021