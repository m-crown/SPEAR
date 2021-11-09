#!/usr/bin/env python
# coding: utf-8

'''
Script to take an input fasta file containing experiments consensus sequences and filters according to %N cutoff.
'''

#requires biopython
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Returns FASTA IDs from a master FASTA file or filtered master FASTA file with %%N content  .')
parser.add_argument('-i', '--input_file', metavar='', type=str,
                    help='Input fasta file with all sequences')
parser.add_argument('-c', '--cutoff', metavar='', type=float, default = 90.5,
                    help='Cutoff value for N content. Default = 90.5')
parser.add_argument('-o', '--output_fasta', metavar='', type=str,
                    help='Destination file for filtered sequences')
parser.add_argument('-l', '--output_list', metavar='', type=str, default = "hoci_seq_ids.txt",
                    help='Destination file for IDs meeting %%N cutoff')

args = parser.parse_args()
print("Running with the following parameters:\nInput file: ", args.input_file,"\nCutoff: ", args.cutoff, "\nOutput FASTA: ", args.output_fasta, "\nOutput_list: ", args.output_list, "\n")

perc_n = {}
low_nsamples = []

short_sequences = []

for record in SeqIO.parse(args.input_file, "fasta"):
    if len(record.seq) == 0:
        short_sequences.append(record)
        perc_n[record.id] = 100
    else:
        seq_length = len(record.seq)
        num_ns = record.seq.count("N")
        perc_nseq = (record.seq.count("N")/len(record.seq)) * 100
        if perc_nseq <= args.cutoff:
            short_sequences.append(record)
            perc_n[record.id] = perc_nseq

print("Found %i sequences with %%N below cutoff" % len(short_sequences))
SeqIO.write(short_sequences, args.output_fasta, "fasta")

with open(args.output_list, 'w') as f:
    for key in perc_n.keys():
        f.write("%s,%s\n" % (key, perc_n[key]))
