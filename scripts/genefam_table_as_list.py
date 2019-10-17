#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys

nffam = sys.argv[1]
nfout = sys.argv[2]
if len(sys.argv)>3:
	a = sys.argv[3]
	if a.isdigit(): writefamname = bool(int(a))
	else: writefamname = bool(eval(a))
else:
	writefamname = False

# could just read input file and write as it reads family blocks
# however non-efficient use of a dict to store family is prefered
# to deal with generic case where families are not sorted in imput file
dfam = {}

with open(nffam, 'r') as ffam:
	for line in ffam:
		fam, prot = line.rstrip('\n').split('\t')
		dfam.setdefault(fam, []).append(prot)

with open(nfout, 'w') as fout:
	for fam, prots in dfam.iteritems():
		if writefamname: l = [fam]+prots
		else: l = prots
		fout.write('\t'.join(l)+'\n')
