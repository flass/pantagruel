#!/usr/bin/python

import sys, os
from ptg_utils import extractCDSFastaFromGFFandGenomicFasta, extractCDSFastaFromGBFF

def checkExistsWithWithoutGz(nf):
	if os.path.exists(nf):
		return nf
	else:
		nfgz = str(nf)+".gz"
	if os.path.exists(nfgz):
		return nfgz
	else:
		return None

nfldirassemb = sys.argv[1]
if len(sys.argv) > 2:
	dirnewcdsfasta = sys.argv[2]
else:
	dirnewcdsfasta = None

with open(nfldirassemb, 'r') as fldirassemb:
	ldirassemb = [line.strip('\n') for line in fldirassemb]


for dirassemb in ldirassemb:
	# parse CDS fasta file
	nfcdsfasta = "%s/%s_cds_from_genomic.fna"%(dirassemb, os.path.basename(dirassemb))
	enfcdsfasta = checkExistsWithWithoutGz(nfcdsfasta)
	if not enfcdsfasta:
		# the '*_cds_from_genomic.fna[.gz]' is missing from the assembly folder (happens for recently published assemblies)
		# will derive its equivalent from the genomic source
		nfgff = "%s/%s_genomic.gff"%(dirassemb, os.path.basename(dirassemb))
		nfgenofna = "%s/%s_genomic.fna"%(dirassemb, os.path.basename(dirassemb))
		nfgbff = "%s/%s_genomic.gbff"%(dirassemb, os.path.basename(dirassemb))
		enfgff =  checkExistsWithWithoutGz(nfgff)
		enfgenofna = checkExistsWithWithoutGz(nfgenofna)
		enfgbff = checkExistsWithWithoutGz(nfgbff)
		if dirnewcdsfasta:
			os.mkdir(os.path.join(dirnewcdsfasta, os.path.basename(dirassemb)))
			nfcdsfasta = "%s/%s/%s_cds_from_genomic.fna.gz"%(dirnewcdsfasta, os.path.basename(dirassemb), os.path.basename(dirassemb))
		if enfgenofna and enfgff:
			# from '*_genomic.gff.gz' and '*_genomic.fna.gz'
			print "create CDS Fasta file extracting data from genomic GFF and Fasta files: '%s' + '%s' -> '%s'"%(nfgff, nfgenofna, nfcdsfasta)
			extractCDSFastaFromGFFandGenomicFasta(enfgff, enfgenofna, nfcdsfasta)
		elif enfgbff:
			# from '*_genomic.gbff.gz'
			print "create CDS Fasta file extracting data from GenBank flat file: '%s' -> '%s'"%(nfgbff, nfcdsfasta)
			extractCDSFastaFromGBFF(enfgbff, nfcdsfasta)
		else:
			missingfiles = '\n'.join(['%s[.gz]'%nf for nf, enf in zip([nfcdsfasta, nfgff, nfgenofna, nfgbff], [enfcdsfasta, enfgff, enfgenofna, enfgbff]) if not enf])
			raise IOError, "could not find any genomic source file to extract CDS from; files are missing:\n"+missingfiles
