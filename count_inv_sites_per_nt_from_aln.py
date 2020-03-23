#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Mar 15 15:32:06 2018
Update to outfile the results Mon Feb 24 10:01:39 2020
@author: gyebra

This script will STDOUT the number of invariable sites for each nucleotide (A, C, G, T).
It will ignore non-ACGT bases so better provide a clean alignment (created e.g. with
snp-sites -c flag).

"""
import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
import pandas as pd

input_file = sys.argv[1]

print("Reading input alignment...")
align = AlignIO.read(input_file, 'fasta') # Whatever format your mutliple sequence alignment is in
summary_align = AlignInfo.SummaryInfo(align)
consensus = summary_align.dumb_consensus()
print("Read alignment of " + str(len(align)) + " sequences and length " + str(align.get_alignment_length()))
print("Calculating PSSM matrix...")
my_pssm = summary_align.pos_specific_score_matrix(consensus)

pssm_list = my_pssm.pssm

dictionary = {}

for i in range(0, len(pssm_list)):
    dictionary[i] = pssm_list[i][1]
print(str(i) + " positions processed...")

pssm_df = pd.DataFrame(dictionary).T

#output_file = input_file + ".tsv"
#pssm_df.to_csv(output_file, sep='\t')

length = float(len(align))
nA = len(pssm_df[pssm_df["A"] == length])
nC = len(pssm_df[pssm_df["C"] == length])
nG = len(pssm_df[pssm_df["G"] == length])
nT = len(pssm_df[pssm_df["T"] == length])

print("Number of invariable sites per nucleotide:")
print("A: " + str(nA))
print("C: " + str(nC))
print("G: " + str(nG))
print("T: " + str(nT))

output_file = input_file + "_invSiteCount.txt"
outfile = open(output_file, 'w')
outfile.write("A\tC\tG\tT\n")
outfile.write(str(nA) + "\t" + str(nC) + "\t" + str(nG) + "\t" + str(nT) + "\n")
outfile.close()


