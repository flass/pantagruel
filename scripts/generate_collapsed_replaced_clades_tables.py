#!/usr/bin/python

import sys, os, glob
from replace_species_by_pop_in_gene_trees import annotatePopulationInSpeciesTree, parseMrBayesConstraints

cladenameprefix = 'clade'
replacedtag = 'leaflabels_Spe2Pop.txt'

dircons = sys.argv[1]
dirrepl = sys.argv[2]
nfoutcons = sys.argv[3]
nfoutrepl = sys.argv[4]

assert os.path.isdir(dircons)
assert os.path.isdir(dirrepl)
foutcons = open(nfoutcons, 'w')
foutrepl = open(nfoutrepl, 'w')

for fam in os.listdir(dircons):
	dirfamcons = os.path.join(dircons, fam)
	lnfcons = [os.path.join(dirfamcons, nf) for nf in os.listdir(dirfamcons)]
	constraintclades = parseMrBayesConstraints(lnfcons)
	for cla, leaflabs in constraintclades.iteritems():
		claid = int(cla.split(cladenameprefix, 1)[-1])
		foutcons.write('\n'.join(["%s\t%d\t%s"%(fam, claid, leaflab) for leaflab in leaflabs])+'\n')

globreplaced = "%s/*%s"%(dirrepl, replacedtag)
lnfreplaced = glob.glob(globreplaced)
for nfreplaced in lnfreplaced:
	fam = os.path.basename(nfreplaced).split('-', 1)[0]
	with open(nfreplaced, 'r') as freplaced:
		for line in freplaced:
			if line.startswith(cladenameprefix): li = line.split(cladenameprefix, 1)[1]
			else: li = line
			foutrepl.write(fam+'\t'+li)
