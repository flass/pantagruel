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
from BCBio import GFF

def parse_strain_info(straininfo_file):
	with open(straininfo_file, 'r') as fin:
		header = fin.readline().rstrip('\n').split('\t')
		dstraininfo = {}
		for line in fin:
			lsp = line.rstrip('\n').split('\t')
			straininfo = dict(zip(header, lsp))
			assid = straininfo['assembly_id']
			dstraininfo[assid] = straininfo
	return dstraininfo

def main(gff_file, fasta_file, straininfo_file=None):
	if fasta_file.endswith('.gz'): ffasta = gzip.open(fasta_file, 'rb')
	else: ffasta = open(fasta_file, 'r')
	fasta_input = SeqIO.to_dict(SeqIO.parse(ffasta, "fasta", alphabet=generic_dna))
	ffasta.close()
	if gff_file.endswith('.gz'):
		fgff = gzip.open(gff_file, 'rb')
		gid = os.path.splitext(os.path.splitext(gff_file)[0])[0]
	else:
		fgff = open(gff_file, 'r')
		gid = os.path.splitext(gff_file)[0]
	genome = list(GFF.parse(fgff, fasta_input))
	fgff.close()
	if straininfo_file: straininfo = parse_strain_info(straininfo_file).get(gid)
	else: straininfo = {}
	# enforce alphabet (bug in GFFParser)
	for seqrec in genome:
		seqrec.seq.alphabet = generic_dna
		if straininfo:
			for feature in seqrec.features:
				if feature.type=='region':
					species = straininfo['genus']+' '+straininfo['species'].strip()
					if not 'organism' in feature.qualifiers:
						sorg = species
						if straininfo['strain']: sorg += ' str. '+straininfo['strain']
						feature.qualifiers['organism'] = [sorg]
					if not 'species' in feature.qualifiers:
						feature.qualifiers['organism'] = [species]
					break # the for feature loop
	# output Genbank
	with open("%s.gbk" % gid, 'w') as gbout_file:
		SeqIO.write(genome, gbout_file, "genbank")
	# output CDS and Proteins sequences
	cdsout_file = open("%s.ffn" % gid, 'w')
	protout_file = open("%s.faa" % gid, 'w')
	for genomerecord in genome:
		for feature in genomerecord.features:
			if feature.type=='CDS':
				qualstr = ''.join(['[%s=%s]'%(qkey, '; '.join(qval)) for qkey, qval in feature.qualifiers.iteritems()])
				product = feature.qualifiers.get('product', '')
				seqcds = feature.location.extract(genomerecord).seq
				cdsout_file.write(">%s %s\n%s\n" % (feature.id, qualstr, str(seqcds)))
				seqprot = seqcds.translate(table=11)
				protout_file.write(">%s %s\n%s\n" % (feature.id, product, str(seqprot)))
	cdsout_file.close()
	protout_file.close()
	# output RNA sequences
	rnaout_file = open("%s.frn" % gid, 'w')
	for genomerecord in genome:
		for feature in genomerecord.features:
			if feature.type.endswith('RNA'):
				qualstr = ' '.join(['[%s=%s]'%(qkey, '; '.join(qval)) for qkey, qval in feature.qualifiers.iteritems()])
				loctag = feature.qualifiers.get('locus_tag', '')
				seqrna = feature.location.extract(genomerecord).seq
				if feature.id is None: feature.id = loctag[0]
				rnaout_file.write(">%s %s\n%s\n" % (feature.id, qualstr, str(seqrna)))
	rnaout_file.close()

if __name__ == "__main__":
	main(*sys.argv[1:])
