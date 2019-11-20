#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, os

gproject = sys.argv[1]
nfgbkin = sys.argv[2]
nfgbkout = sys.argv[3]
nfstraininfo = sys.argv[4]

fgbkin = open(nfgbkin, 'r')
fgbkout = open(nfgbkout, 'w')

dstraininfo = {}
with open(nfstraininfo, 'r') as fstraininfo:
	header = fstraininfo.readline().rstrip('\n').split('\t')
	gprojfield = header.index('assembly_id')
	for line in fstraininfo:
		lsp = line.rstrip('\n').split('\t')
		gproj = lsp[gprojfield]
		dstraininfo[gproj] = dict(zip(header, lsp))

straininfo = dstraininfo[gproject]

for line in fgbkin:
	fgbkout.write(line)
	if line.startswith('                     /strain'):
		strain = line.rstrip(' \n').split('=')[1].strip('"').replace(' ', '')
		#~ try:
			#~ straininfo = dstraininfo[strain]
		#~ except KeyError, e:
			#~ raise KeyError, "missing information on strain '%s'; please complete info in file '%s'"%(strain, nfstraininfo)
		# add a line with the taxid
		taxidline =    '                     /db_xref="taxon:%s"\n'%straininfo['taxid']
		fgbkout.write(taxidline)
		
fgbkin.close()
fgbkout.close()
