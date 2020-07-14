#!/usr/local/env python

"""
This script strips regions of a fasta alignment based on the 
provided coordinates 

Usage:
	python subsetAln.py input_alignment output_alignment coordinates

G. Yebra, created on August 2017
"""

# import packages
from Bio import AlignIO
import sys

# parse positional arguments
infile = sys.argv[1]
outfile = sys.argv[2]
coords = sys.argv[3]


alignment = AlignIO.read(infile, "fasta")
coordinates = open(coords, "r")
lines=coordinates.readlines()
start = 0
sliced_aln = alignment[:,:0]
for line in lines:
    end=int(line.split('\t')[0])
    print(end)
    sliced_aln = sliced_aln + alignment[:,start:end] 
    start=int(line.split('\t')[1])-1
    print(start)
sliced_aln += alignment[:,start:] 
AlignIO.write(sliced_aln, outfile,"fasta")
coordinates.close()

