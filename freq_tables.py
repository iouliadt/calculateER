import sys
import time
from Bio import SeqIO

start_time = time.time()
# source files
gbPath = input("Type genbank file path\n")
fastaPath = input("Type fasta file path\n")
# parse the genbank file
gbFile = SeqIO.parse(gbPath, "genbank")
# parse fasta file
fastaFile = SeqIO.parse(fastaPath, "fasta")

# Extract gene coords & gene names from genbank file
geneLocation = []
geneName = []
for record in gbFile:
    for f in record.features:
        if f.type == "gene":
            geneLocation.append(f.location)
            if "gene" in f.qualifiers:
                geneName.append(f.qualifiers["gene"])

# ! create empty array to attach to the following dictionary
empty = []
for i in range(len(geneName)):
    # ! Flatten the data tree
    geneName[i] = geneName[i][0]
    # ! Workaround for the creation of the dictionary
    empty.append([])
# ! create dict with gene name and empty lists (in order to append the coords)
genes = dict(zip(geneName, empty))

# slice fasta according to the extracted coords and store to the dictionary "genes"
for seq_record in fastaFile:
    for i in range(len(geneLocation)):
        genes[geneName[i]].append(
            seq_record.seq[geneLocation[i].start: geneLocation[i].end]
        )

# Create a list from the dictionary keys (the gene names)
names = list(genes.keys())
# ! print(genes[geneName[5]][5][0]) PRINTS A SINGLE NUCLEOTIDE


# Create lists of column values of the genes
def nucLists(gene):
    nucCols = []
    for i in range(len(gene[0])):
        x = ""
        for j in range(1, len(gene)):
            x = x + gene[j][i]
        nucCols.append(x)
    return nucCols


# count gaps, A, T, C, Gs and write header to file
def nucTable_gaps(nucCols_ORF, ORF, ORF_name):
    with open(
        "{}_nucfreq_gaps.txt".format(
            ORF_name
        ),
        "a",
    ) as f:
        f.write("NucSite\tRefNuc\tRefGaps\tAcnt\tTcnt\tGcnt\tCcnt\n")
        for i in range(len(nucCols_ORF)):
            print(
                i + 1,
                "\t",
                ORF[0][i],
                "\t",
                ORF[0][i].count("-"),
                "\t",
                nucCols_ORF[i].count("A"),
                "\t",
                nucCols_ORF[i].count("T"),
                "\t",
                nucCols_ORF[i].count("G"),
                "\t",
                nucCols_ORF[i].count("C"),
                file=f,
            )


# ask for user input to create a file per ORF requested
for i in range(len(names)):
    print("For", names[i], "type", i)
x = int(input())

# Create the column lists
nucCols = nucLists(genes[geneName[x]])
# Write the tables
nucTable_gaps(nucCols, genes[geneName[x]], names[x])


# Custom translate function (to replace "---" with "-")
def translate_dna(sequence):
    codontable = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGC": "C",
        "TGT": "C",
        "TGA": "*",
        "TGG": "W",
        "---": "-",
    }

    prot = []
    for n in range(0, len(sequence), 3):
        if sequence[n: n + 3] in codontable:
            residue = codontable[sequence[n: n + 3]]
        else:
            residue = "X"
        prot.append(residue)
    return "".join(prot)


translations = []
for i in range(len(genes[geneName[x]][0])):
    translations.append(translate_dna(genes[geneName[x]][i]))

polarity = {
    "T": "polar",
    "S": "polar",
    "N": "polar",
    "Q": "polar",
    "Y": "polar",
    "K": "basic polar",
    "R": "basic polar",
    "H": "basic polar",
    "D": "acidic polar",
    "E": "acidic polar",
    "A": "nonpolar",
    "L": "nonpolar",
    "I": "nonpolar",
    "V": "nonpolar",
    "G": "nonpolar",
    "P": "nonpolar",
    "W": "nonpolar",
    "M": "nonpolar",
    "C": "nonpolar",
    "F": "nonpolar",
    "*": "",
}


def aaLists(translatedGene):
    # aaCols = []
    ref_aa = []
    alt_aa = []
    site = []
    mutation = []
    ref_polarity = []
    alt_polarity = []
    for i in range(len(translatedGene[0])):
        x = 0
        for j in range(1, len(translatedGene)):
            if translatedGene[0][i] == translatedGene[j][i]:
                continue
            else:
                x = len(ref_aa)
                if translatedGene[0][i] != "X" and translatedGene[j][i] != "X":
                    ref_aa.append(translatedGene[0][i])
                    alt_aa.append(translatedGene[j][i])
                    site.append(i + 1)
                    mut_string = ref_aa[x] + str(site[x]) + alt_aa[x]
                    mutation.append(mut_string)
                    if ref_aa[x] != "*":  # ref_aa[x] != 'X':
                        refpol = polarity[ref_aa[x]]
                        ref_polarity.append(refpol)
                    else:
                        refpol_stop = "stop codon"
                        ref_polarity.append(refpol_stop)
                    if alt_aa[x] != "*" and alt_aa[x] != "-":
                        altpol = polarity[alt_aa[x]]
                        alt_polarity.append(altpol)
                    else:
                        altpol_stop = "stop codon"
                        alt_polarity.append(altpol_stop)

    dict = {
        "Mutation": mutation,
        "Ref_aa": ref_aa,
        "Ref_polarity": ref_polarity,
        "aa_position": site,
        "Query_aa": alt_aa,
        "Alt_polarity": alt_polarity,
    }
    # Print the dictionary
    counter = 0
    for keys in dict:
        if counter != len(dict)-1:
            print(keys, end=" ")
        else:
            print(keys)
        counter += 1
    for i in range(len(dict["Mutation"])):
        counter = 0
        for j in dict:
            if counter != len(dict)-1:
                print(dict[j][i], end=" ")
            else:
                print(dict[j][i])
            counter += 1


# print(translations)
aaLists(translations)

print("\n--- %s seconds ---" % (time.time() - start_time))
sys.exit()