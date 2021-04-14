#!/usr/bin/python2.7
"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
GFFandFasta2GenBank.py <GFF annotation file> <FASTA sequence file>

Credit Brad Chapman https://www.biostars.org/p/2492/
"""

import sys
import os

from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from BCBio import GFF

def main(gff_file, fasta_file):
    out_file = "%s.gbk" % os.path.splitext(gff_file)[0]
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", DNAAlphabet)
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(gff_iter, out_file, "genbank")

if __name__ == "__main__":
    main(*sys.argv[1:])
