# python calculateER.py -f multifasta_example.fasta -d1 23/3/2021 -d2 6/5/2021

from Bio import SeqIO
import argparse
from datetime import date, datetime
import sys

parser = argparse.ArgumentParser(
        description='Calculate the rate of evolution between two time points.')

parser.add_argument('-f', '--fasta', required=True,
                    help='Path to multifasta file with the two consensus sequences')
parser.add_argument('-d1', '--date1', required=True,
                    help='Respective date of the first sequence in the format of: dd/mm/yyyy')
parser.add_argument('-d2', '--date2', required=True,
                    help='Respective date of the second sequence in the format of: dd/mm/yyyy')

# Display help message when svartools is called with wrong arguments
args = parser.parse_args()

if not vars(args) or len(sys.argv) <= 2:
    parser.print_help()
    parser.exit(1)

# Parse fasta file
records = list(SeqIO.parse(args.fasta, "fasta"))
seq1 = list(records[0].seq)
seq2 = list(records[1].seq)

# Check if the sequences have the same length
if len(seq1) != len(seq2):
    raise ValueError("Sequences of different length!")

total_sites = len(seq1)
changes = 0

# Count the nucleotide changes between the two sequences
for nuc in range(len(seq1)):
    if seq1[nuc] != seq2[nuc]:
        if seq1[nuc] != "N" and seq2[nuc] != "N" and seq1[nuc] != "-" and seq2[nuc] != "-":
            changes += 1

print("The number of nucleotide changes between the two sequences, ignoring gaps and N characters, is {}".format(changes))

# Calculate the substitutions per site
subPerSite = changes / total_sites

date1 = args.date1
date2 = args.date2

date_format = '%d/%m/%Y'

try:
    d1 = datetime.strptime(args.date1, date_format)
    d2 = datetime.strptime(args.date2, date_format)
except ValueError:
    raise ValueError("Incorrect data format, should be dd/mm/yyyy")


# Calculate the number of days between the two dates
delta = d2 - d1
days = delta.days

# Divide sub per site by the days to get sub per site per day
subPerSitePerday = subPerSite/days
# Multiply by 365 to get subs per site per year
subPerSitePerYear = subPerSitePerday * 365

print("The average number of substitutions per site per year is {}".format(subPerSitePerYear))
