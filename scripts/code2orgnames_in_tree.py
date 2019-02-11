#!/usr/bin/python

import tree2
import sys

sys.setrecursionlimit(20000)

nftreein = sys.argv[1]
nfref = sys.argv[2]
nfout = sys.argv[3]

t = tree2.read_check_newick(nftreein, maxNoBS=-1)

with open(nfref, 'r') as fref:
	for line in fref:
		lsp = line.rstrip('\n').split('\t')
		code = lsp[0]
		orga = lsp[1]
		stra = lsp[2]
		if stra:
			stra = stra.split('type strain:')[-1]
			if ' = ' in stra: stra = stra.split(' = ')[0]
			if '; ' in stra: stra = stra.split('; ')[0]
		if ' = ' in orga: orga = orga.split(' = ')[0]
		if '; ' in orga: orga = orga.split('; ')[0]
		
		#~ if orga.endswith(' sp'): org +='.'
		if not (stra in orga):
			orga += ' str. '+stra
		elif stra:
			if ' str. ' not in orga:
				orga = orga.replace(stra, 'str. '+stra)
		
		orga = orga.replace('(', '').replace(')', '').replace('.', '').replace(':', '_').replace(' ', '_')
		orga += '__'+code
		print code, orga
		node = t[code]
		if node:
			t[code].edit_label(orga)
			print orga
		else:
			print "could not find %s in tree"%code

t.write_newick(nfout)
