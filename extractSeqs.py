#!/usr/bin/env python

"""
This script extract sequences from an alignment whose name is in a list 
of names provided as text file

Usage: 
	python extractSeqs.py list.txt aln.fa

G. Yebra, created on August 2017
"""

import sys
sample_list = sys.argv[1]
aln = sys.argv[2]

from Bio import AlignIO
alignment = AlignIO.read(aln, "fasta")
unq_list = open(sample_list, "r")
unq_seqs = unq_list.readlines()

outfile_name = aln + ".filtered"
outfile = open(outfile_name, 'w')

for unq_seq in unq_seqs:
    for record in alignment:
        if record.id == unq_seq.replace("\n",""):
            outfile.write(">" + record.id + "\n" + str(record.seq) + "\n") 
outfile.close
unq_list.close
            