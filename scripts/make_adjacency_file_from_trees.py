#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import tree2
import sys

nflnfgt = sys.argv[1]
nfout = sys.argv[2]

if len(sys.argv) > 3:
	maxggdist = int(sys.argv[3])
else:
	maxggdist = 1
with open(nflnfgt,'r') as flnfgt:
	lnfgt = [line.rstrip('\n') for line in flnfgt]

# store gene labels
dspe2genes = {}
ngenes = 0
for nfgt in lnfgt:
	try:
		gt = tree2.read_newick(nfgt)
	except ValueError:
		# assume bootstrap/bayesian samples don't have branch lengths
		gt = tree2.read_multiple_newick(nfgt, branch_lengths=False)[0]
	genes = gt.get_leaf_labels()
	for g in genes:
		# assume gene labels to follow the syntax SPECIES_1234
		s, n = g.split('_')
		dspe2genes.setdefault(s, []).append(int(n))
	ngenes += len(genes)

lspe = dspe2genes.keys()
lspe.sort()
print "%d species: %s, ...; %d gene families; %d genes."%(len(lspe), ', '.join(lspe[:5]), len(lnfgt), ngenes)

# make adjacency file 
with open(nfout, 'w') as fout:
	for s in lspe:
		lg = dspe2genes[s]
		lg.sort()
		for k in range(len(lg)-1):
			a = lg[k]
			b = lg[k+1]
			if (b - a) <= maxggdist:
				# only consider single-replicon, linear topology genome
				# have to query db to take into account multi-replicon and possible circular topology
				fout.write("%s_%d %s_%d\n"%(s, a, s, b))

