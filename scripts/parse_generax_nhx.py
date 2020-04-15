#!/usr/bin/env python2.7

import tree2

import sys, os, getopt

#sys.setrecursionlimit(20000)

#deventshort = {'S':'S', 'D':'D', 'H':'T'}

def parse_branch_annot(annot):
	levt = []
	for elt in annot.split(':'):
		if not elt=='&&NHX':
			evtype, evloc = elt.split('=')
			if evloc!='N':
				if (elt[0] in ('S', 'D')):
					levt.append((evtype, evloc, ''))
				elif elt[0]=='H':
					don, rec = evloc.split('@')[1:]
					levt.append(('T', rec, don))
	return levt
	

nfinnhx = sys.argv[1]
nfout = sys.argv[2]

lgt = read_multiple_newick(nfinnhx, keep_comments=True)

devtcount = {}
for gt in lgt:
	for node in gt:
		for evt in parse_branch_annot(node.comment()):
			devtcount[evt] = devtcount.setdefault(evt, 0) + 1

with open(nfout, 'w') as fout:
	for evt, n in devtcount.iteritems():
		fout.write('\t'.join(evt+(str(n,))+'\n')
