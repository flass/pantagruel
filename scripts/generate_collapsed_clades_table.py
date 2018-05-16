#!/usr/bin/python

import sys, os
from replace_species_by_pop_in_gene_trees import annotatePopulationInSpeciesTree, parseMrBayesConstraints

cladenameprefix = 'clade'

dircons = sys.argv[1]
nfout = sys.argv[2]

fout = open(nfout, 'w')
for fam in os.listdir(dircons):
	dirfamcons = os.path.join(dircons, fam)
	lnfcons = [os.path.join(dirfamcons, nf) for nf in os.listdir(dirfamcons)]
	constraintclades = parseMrBayesConstraints(lnfcons)
	for cla, leaflabs in constraintclades.iteritems():
		claid = int(cla.split(cladenameprefix, 1)[-1])
		fout.write('\n'.join(["%s\t%d\t%s"%(fam, claid, leaflab) for leaflab in leaflabs])+'\n')
			
