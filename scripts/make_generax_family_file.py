#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, os, glob, getopt
import tree2

def usage():
	s = '[need help message here]'
	return s

opts, args = getopt.getopt(sys.argv[1:], 'a:o:hv', ['alignments=', 'gene-trees=', 'model=', 'out='])
dopt = dict(opts)
if ('-h' in dopt) or ('--help' in dopt):
	print usage()
	sys.exit(0)

gttag = '-Gtree.nwk'
alitag = 'aln'

dirali = dopt.get('-a', dopt.get('--alignments'))
dirgt = dopt.get('--gene-trees')
model = dopt.get('model', 'GTR+G')
nfout = dopt.get('-o', dopt.get('--out'))

if not nfout:
	raise ValueError, "you must specify an ouput file with -o|--out"
if not dirali:
	raise ValueError, "you must specify an input alignment folder with -a|--alignments"

lnfali = glob.glob('%s/*%s'%(dirali, gttag))
lnfgt = glob.glob('%s/*%s'%(dirgt, gttag))

fout = open(nfout, 'w')
fout.write('[FAMILIES]\n')

for nfali in lnfali:
	bnali = os.path.basename(nfali)
#	dirgt = os.path.dirname(nfgt)
	fam = bnali.rsplit('.')[0]
	if dirgt:
		nfgt = glob.glob("%s/*%s*%s"%(dirgt, fam, gttag))
		bnmap = bnrad2+'.link'
		nfmap = os.path.join(dirgt, bnmap)
		with open(nfmap, 'w') as fmap:
			# generate the mapping file
			gt = tree2.read_newick(nfgt)
			for leaflab in gt.get_leaf_labels():
				fmap.write('%s\t%s\n'%(leaflab, leaflab.split('_', 1)[0]))
	else:
		nfmap = None
		nfgt = None
	
	# write out
	fout.write('- %s\n'%fam)
	fout.write('alignment = %s\n'%(os.path.join(dirgt, bnaln)))
	fout.write('subst_model = %s\n'%(model))
	if nfgt: fout.write('starting_gene_tree = %s\n'%(nfgt))
	if nfmap: fout.write('mapping = %s\n'%(os.path.join(dirgt, bnmap)))
	
fout.close()