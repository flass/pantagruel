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


def main(nflnfgbff, nfrnaout):
	frnaout = open(nfrnaout, 'w')
	with open(nflnfgbff, 'r') as flnfgbff:
		lnfgbff = [line.strip('\n') for line in flnfgbff]
	for nfgbff in lnfgbff:
		if nfgbff.endswith('.gz'):
			fgbff = gzip.open(nfgbff, 'rb')
		else:
			fgbff = open(nfgbff, 'r')
		genome = list(SeqIO.parse(fgbff, format="genbank"))
		fgbff.close()
		organism = None
		# output RNA sequences
		for genomerecord in genome:
			grid = genomerecord.id
			for feature in genomerecord.features:
				if (feature.type == 'source') and (organism is None):
					organism = feature.qualifiers.get('organism', [''])[0]
					strain = feature.qualifiers.get('strain', [''])[0]
					print 'grid, organism, strain:', grid, organism, strain
				if feature.type == 'rRNA':
					product = feature.qualifiers.get('product', [''])
					if '16S ribosomal RNA' in product:
						loctag = feature.qualifiers.get('locus_tag', '')
						seqrna = feature.location.extract(genomerecord).seq
						if len(seqrna)<1000:
							continue
						qualstr = ''.join(['[%s=%s]'%(qkey, '; '.join(qval)) for qkey, qval in feature.qualifiers.iteritems()])
						idstr = (organism+('_'+strain if len(organism.split())<=2 else '')+'_'+grid).replace(' ', '_')
						print idstr
						frnaout.write(">%s %s\n%s\n" % (idstr, qualstr, str(seqrna)))
						break
			else:
				continue
			break
	frnaout.close()

if __name__ == "__main__":
	main(*sys.argv[1:])
