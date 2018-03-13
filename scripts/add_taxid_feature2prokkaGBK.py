#!/usr/bin/python

import sys, os

nfgbkin = sys.argv[1]
nfgbkout = sys.argv[2]

fgffin = open(nfgffin, 'r')
fgffout = open(nfgffout, 'w')

with open(nfstraininfo, 'r') as fstraininfo:
	header = fstraininfo.readline().rstrip('\n').split('\t')
	for line in fstraininfo:
		lsp = line.rstrip('\n').split('\t')
		strain = lsp[3]
		dstraininfo[strain] = dict(zip(header, lsp))

for line in fgbkin:
	fgbkout.write(line)
	if line.startswith('                     /strain'):
		strain = line.rstrip(' \n').split('=')[1].strip('"')
		straininfo = dstraininfo[strain]
		# add a line with the taxid
		taxidline =    '                     /db_xref="taxon:%s"\n'%straininfo['taxid']
		fgbkout.write(taxidline)
		
fgbkin.close()
fgbkout.close()
