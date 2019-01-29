#!/usr/bin/python2.7
"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
GFFGenomeFasta2GenBankCDSProtFasta.py <GFF annotation file> <FASTA sequence file>

Credit Brad Chapman https://www.biostars.org/p/2492/
"""

import sys
import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from BCBio import GFF

def main(gff_file, fasta_file):
	fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", alphabet=generic_dna))
	genome = list(GFF.parse(gff_file, fasta_input))
	# output Genbank
	with open("%s.gbk" % os.path.splitext(gff_file)[0], 'w') as gbout_file:
		SeqIO.write(genome, gbout_file, "genbank")
	# output CDS and Proteins
	cdsout_file = open("%s.ffn" % os.path.splitext(gff_file)[0], 'w')
	protout_file = open("%s.faa" % os.path.splitext(gff_file)[0], 'w')
	for genomerecord in genome:
		for feature in genomerecord.features:
			if feature.type=='CDS':
				product = feature.qualifiers.get('product', '')
				seqcds = feature.location.extract(genomerecord).seq
				cdsout_file.write(">% %s\n%s\n" % (feature.id, product, str(seqcds)))
				seqprot = seqcds.translate(table=11)
				protout_file.write(">% %s\n%s\n" % (feature.id, product, str(seqprot)))
	cdsout_file.close()
	protout_file.close()

if __name__ == "__main__":
	main(*sys.argv[1:])
