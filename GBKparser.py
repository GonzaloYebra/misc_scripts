#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script takes a list of locus tags and a genbank file and will extract
two fasta files, one containing the nucleotide sequences and one the amino 
acid sequences corresponding of the genes identified with those locus tags.
It will also output a csv file with information about the genes.

Usage:
    python GBKparser.py [options] locus_tag_list genbank_file.gbk

Author: G. Yebra
Created on 06/12/2018
"""

# import packages
import argparse
import pandas as pd
from Bio import SeqIO

# parse arguments and flags
parser = argparse.ArgumentParser()

parser.add_argument('tags',type=str,help='locus tag list file')
parser.add_argument('genbank',type=str,help='genbank file')
parser.add_argument("-nt","--nucleotides", help="extract nucleotide sequences", action="store_true")
parser.add_argument("-aa","--protein", help="extract amino acid sequences", action="store_true")
parser.add_argument("-f","--features", help="extract gene features in a csv file", action="store_true")

args = parser.parse_args()

# create output files
if args.nucleotides:
    outfile_nt = open("nt_seqs.fna", 'w') 
if args.protein:
    outfile_aa = open("aa_seqs.faa", 'w') 

# loop over the list of locus_tag and extract the nt and aa sequences and 
# their data from the genbank file
dct = {}
for tag in open(args.tags):
    for record in SeqIO.parse(args.genbank, "genbank"):
        for f in record.features:
            if f.type == "CDS" and "locus_tag" in f.qualifiers:
                locus_tag = f.qualifiers["locus_tag"][0]
                if locus_tag  == tag.strip('\n'):
                    # extract nt sequence in fasta format
                    if args.nucleotides:
                        outfile_nt.write(">" + f.qualifiers["locus_tag"][0] + "\n" + str(f.location.extract(record.seq)) + "\n")
                    # extract aa sequence in fasta format
                    if args.protein:
                        outfile_aa.write(">" + f.qualifiers["locus_tag"][0] + "\n" + f.qualifiers["translation"][0] + "\n") 
                    # store the genes' features
                    if args.features:
                        locus_product = f.qualifiers["product"][0]
                        locus_start = str(f.location.start)
                        locus_end = str(f.location.end)
                        locus_strand = str(f.location.strand)
                        if "gene" in f.qualifiers:
                            gene_label = f.qualifiers["gene"][0]
                        else:
                            gene_label = "NA"
                            
    # add to the dictionary the gene's characteristics
    if args.features:
        dct[tag.strip('\n')] = {"prot_start": locus_start, "prot_end": locus_end, "prot_product": locus_product, "prot_len": int(locus_end) - int(locus_start), "prot_strand": locus_strand, "gene": gene_label}

# close sequence files    
if args.nucleotides:
    outfile_nt.close()
if args.protein:
    outfile_aa.close()

# create dataframe from dictionary and export as csv file
if args.features:
    df = pd.DataFrame.from_dict(dct).T
    df.to_csv("features.csv", sep=',')
