#!/usr/bin/env python3

from Bio import SeqIO
from sys import argv

input_fasta = argv[1]

if len(argv) != 2:
    print("Usage: python3 split_fastas.py <input_fasta>")
    exit()


all_sequences = list(SeqIO.parse(input_fasta, "fasta"))

non_dated_sequences = [seq for seq in all_sequences if seq.id.endswith("_ND")]
dated_sequences = [seq for seq in all_sequences if not seq.id.endswith("_ND")]

for nd_seq in non_dated_sequences:
    output_file = f"{nd_seq.id}.fasta"
    with open (output_file, "w") as f:
        # write the non-dated sequences to the file
        SeqIO.write(nd_seq, f, "fasta")
        # write the dated sequences to the file
        SeqIO.write(dated_sequences, f, "fasta")





