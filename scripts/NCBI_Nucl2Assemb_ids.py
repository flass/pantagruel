#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, os
from Bio import Entrez

nfnucid = sys.argv[1]
Entrez.email =  sys.argv[2]

with open(nfnucid, 'r') as fnucid:
	lnucid = [line.rstrip('\n') for line in fnucid]
	#~ header = fnucid.readline().rstrip('\n').split('\t')
	#~ ldnucid = [dict(zip(header, line.rstrip('\n').split('\t'))) for line in fnucid]

for nucid in lnucid:
#~ for dnucid in ldnucid:
	#~ nucid = dnucid['Nuc_accession_long']
	print nucid,
	#~ if nucid.startswith('NZ_'):
		#~ nucid = nucid.split('_', 1)[0]
	try:
		erecord = Entrez.read(Entrez.elink(db="assembly", dbfrom="nucleotide", id=nucid))
		assUid = erecord[0]['LinkSetDb'][0]['Link'][0]['Id']
		esum = Entrez.read(Entrez.esummary(db='assembly', id=assUid))
		print esum['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
	except IndexError, e:
		#~ bioprojid = dnucid['Bioproject']
		#~ erecord = Entrez.read(Entrez.elink(db="assembly", dbfrom="bioproject", id=bioprojid, idtype='acc'))
		#~ assUid = erecord[0]['LinkSetDb'][0]['Link'][0]['Id']
		#~ esum = Entrez.read(Entrez.esummary(db='assembly', id=assUid))
		#~ print esum['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
		print ''
