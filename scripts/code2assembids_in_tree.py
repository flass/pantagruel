#!/usr/bin/python

import tree2
import sys

nftreein = sys.argv[1]
nfref = sys.argv[2]
nfout = sys.argv[3]

t = tree2.Node(file=nftreein, unrooted=True)

with open(nfref, 'r') as fref:
	for line in fref:
		lsp = line.rstrip('\n').split('\t')
		ass = lsp[0]
		code = lsp[1]
		ass += '__'+code
		print code, ass
		node = t[code]
		if node:
			t[code].edit_label(ass)
			print ass
		else:
			print "could not find %s in tree"%code

t.write_newick(nfout)
