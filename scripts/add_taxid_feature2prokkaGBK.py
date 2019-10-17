#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, os

nfgbkin = sys.argv[1]
nfgbkout = sys.argv[2]
nfstraininfo = sys.argv[3]

fgbkin = open(nfgbkin, 'r')
fgbkout = open(nfgbkout, 'w')

dstraininfo = {}
with open(nfstraininfo, 'r') as fstraininfo:
	header = fstraininfo.readline().rstrip('\n').split('\t')
	strainfield = header.index('strain')
	for line in fstraininfo:
		lsp = line.rstrip('\n').split('\t')
		strain = lsp[strainfield]
		dstraininfo[strain] = dict(zip(header, lsp))

for line in fgbkin:
	fgbkout.write(line)
	if line.startswith('                     /strain'):
		strain = line.rstrip(' \n').split('=')[1].strip('"').replace(' ', '')
		try:
			straininfo = dstraininfo[strain]
		except IndexError, e:
			raise IndexError, "missing information on strain '%s'; please complete info in file '%s'"%(strain, nfstraininfo)
		# add a line with the taxid
		taxidline =    '                     /db_xref="taxon:%s"\n'%straininfo['taxid']
		fgbkout.write(taxidline)
		
fgbkin.close()
fgbkout.close()
