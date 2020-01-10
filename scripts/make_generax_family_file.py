#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, os, glob
import tree2

gttag = '-Gtree.nwk'

dirgt = sys.argv[1]
nfout = sys.argv[2]
if len(sys.argv) > 3:
	model = sys.argv[3]
else:
	model = 'GTR+G'

#with open(nflnfgt, 'r') as flnfgt:
#	lnfgt = [line.rstrip('\n') for line in flnfgt]
lnfgt = glob.glob('%s/*%s'%(dirgt, gttag))

fout = open(nfout, 'w')
fout.write('[FAMILIES]\n')

for nfgt in lnfgt:
	bngt = os.path.basename(nfgt)
#	dirgt = os.path.dirname(nfgt)
	bnrad1 = bngt.rsplit('.', 1)[0]
	bnrad2 = bnrad1.rsplit('-', 1)[0]
	if bnrad2.endswith('-replaced'): bnrad3 = bnrad2.rsplit('-', 1)[0]
	else: bnrad3 = bnrad2
	fam = bnrad1.split('-')[0]
	# assume everything is in the same folder
	bnaln = bnrad2+'.aln'
#	bnst = bnrad3+'-Stree.nwk'
	bnmap = bnrad2+'.link'
	nfmap = os.path.join(dirgt, bnmap)
	with open(nfmap, 'w'):
		# generate the mapping file
		gt = tree2.read_newick(nfgt)
		for leaflab in gt.get_leaf_labels():
			nfmap.write('%s\t%s\n'%(leaflab, leaflab.split('_', 1)[0]))
	
	# write out
	fout.write('- %s\n'%fam)
	fout.write('mapping = %s\n'%(os.path.join(dirgt, bnmap)))
	fout.write('alignment = %s\n'%(os.path.join(dirgt, bnaln)))
	fout.write('subst_model = %s\n'%(model))
	fout.write('starting_gene_tree = %s\n'%(nfgt))
	
fout.close()