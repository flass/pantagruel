#!/usr/bin/python

import sys

nffam = sys.argv[1]
nfout = sys.argv[2]
writefamname = bool(sys.argv[3])

dfam = {}

with open(nffam, 'r') as ffam:
	fam, prot = line.rstrip('\n').split('\t')
	dfam.setdefault(fam, []).append(prot)

with open(nfout, 'w') as fout:
	for fam, prots in dfam.iteritems():
		if writefamname: l = [fam]+prots
		else: l = prots
		fout.write('\t'.join(l)+'\n')
