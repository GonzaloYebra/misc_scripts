#!/usr/local/env python

"""
This script subsets a region of a fasta alignment based on the provided coordinates 

Usage:
	python subsetAln.py inputfile outputfile startposition endposition

G. Yebra, created on August 2017
"""

# import packages
from Bio import AlignIO
from Bio import SeqIO
import sys

# parse positional arguments
infile = sys.argv[1]
outfile = sys.argv[2]
start = int(sys.argv[3])
stop = int(sys.argv[4])

# subset alignment
align = AlignIO.read(infile, "fasta")
records = align[:,start+1:stop+1]

# export output alignment
SeqIO.write(records, outfileName, "fasta")
