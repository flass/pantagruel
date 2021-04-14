#!/usr/bin/python2.7
"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
GFFGenomeFasta2GenBankCDSProtFasta.py <GFF annotation file> <FASTA sequence file>

Credit Brad Chapman https://www.biostars.org/p/2492/
"""

import sys
import os
import gzip

from Bio import SeqIO
from Bio.Alphabet import generic_dna


def main(gbff_file):
	if gbff_file.endswith('.gz'):
		fgbff = gzip.open(gbff_file, 'rb')
		gid = os.path.splitext(os.path.splitext(gbff_file)[0])[0]
	else:
		fgbff = open(gbff_file, 'r')
		gid = os.path.splitext(gbff_file)[0]
	genome = list(SeqIO.parse(fgbff, format="genbank"))
	fgbff.close()
	# output RNA sequences
	rnaout_file = open("%s.frn" % gid, 'w')
	for genomerecord in genome:
		for feature in genomerecord.features:
			if feature.type.endswith('RNA'):
				product = feature.qualifiers.get('product', '')
				loctag = feature.qualifiers.get('locus_tag', '')
				seqrna = feature.location.extract(genomerecord).seq
				if feature.id in [None, '<unknown id>']:
					feature.id = loctag[0]
				qualstr = ''.join(['[%s=%s]'%(qkey, '; '.join(qval)) for qkey, qval in feature.qualifiers.iteritems()])
				rnaout_file.write(">%s %s\n%s\n" % (feature.id, qualstr, str(seqrna)))
	rnaout_file.close()

if __name__ == "__main__":
	main(sys.argv[1])
