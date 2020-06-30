#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys, os, glob, getopt
import tree2

print "# call was:", ' '.join(sys.argv)

def usage():
	s = '[need help message here]'
	return s

opts, args = getopt.getopt(sys.argv[1:], 'l:a:o:hv', ['alignment_list=', 'alignments=', 'gene-trees=', 'out=', \
													'model=', 'per-family', 'skip-abs-gt', \
													'alitag=', 'gttag=', 'sttag=', 'gftag=', \
													'help', 'verbose'])
dopt = dict(opts)
if ('-h' in dopt) or ('--help' in dopt):
	print usage()
	sys.exit(0)

verbose = (('-v' in dopt) or ('--verbose' in dopt))

gttag = dopt.get('--gttag', '-Gtree.nwk')
sttag = dopt.get('--sttag', '-Stree.nwk')
gftag = dopt.get('--gftag', '.generax_families')
alitag = dopt.get('--alitag', '.aln')
skipabsgt = ('--skip-abs-gt' in dopt)

dirali = dopt.get('-a', dopt.get('--alignments'))
nflistali = dopt.get('-l', dopt.get('--alignment_list'))
if not (dirali or nflistali):
	raise ValueError, "you must specify an input alignment list with -l|--alignment_list and/or an alignment folder with -a|--alignments"
dirgt = dopt.get('--gene-trees')
model = dopt.get('--model', 'GTR+G')
perfam = ('--per-family' in dopt)
nfout = dopt.get('-o', dopt.get('--out'))
if not nfout:
	raise ValueError, "you must specify an ouput file/folder with -o|--out"
if not perfam:
	fout = open(nfout, 'w')
	fout.write('[FAMILIES]\n')
else:
	dirout = nfout
	if not os.path.isdir(dirout):
		raise ValueError, "specified output directory '%s' does not exist / is not a directory / cannot be accessed"%dirout

lnfali = []
if dirali:
	lnfali += glob.glob('%s/*%s'%(dirali, alitag))
if nflistali:
	with open(nflistali) as flistali:
		lnfali += [line.rstrip('\n') for line in flistali]

if dirgt and skipabsgt:
	lnfgt = glob.glob('%s/*%s'%(dirgt, gttag))
	sgtfam = set([os.path.basename(nfgt).split('.')[0].split('-')[0] for nfgt in lnfgt])

for nfali in lnfali:
	bnali = os.path.basename(nfali)
	fam = bnali.split('.')[0].split('-')[0]
	bnrad = bnali.rsplit('.',1)[0]
	if dirgt:
		if skipabsgt:
			if not (fam in sgtfam):
				if verbose: print "skip gene family %s in the absence of a matching gene tree in folder %s"%(fam, dirgt)
				continue # the for nfali loop	
		gfgt = "%s/*%s*%s"%(dirgt, bnrad, gttag)
		lnfgt = glob.glob(gfgt)
		if not lnfgt:
			if skipabsgt:
				if verbose: print "no file match the pattern: %s; skip gene family %s in the absence of a gene tree"%(gfgt, fam)
				continue # the for nfali loop
			else:
				raise IndexError, "found no gene tree file matching the pattern: %s"%gfgt
		nfgt = lnfgt[0]
		bnmap = bnrad+'.link'
		nfmap = os.path.join(dirgt, bnmap)
		with open(nfmap, 'w') as fmap:
			# generate the mapping file
			gt = tree2.read_newick(nfgt)
			for leaflab in gt.get_leaf_labels():
				fmap.write('%s\t%s\n'%(leaflab, leaflab.split('_', 1)[0]))
	else:
		nfmap = None
		nfgt = None
	
	if perfam:
		fout = open(os.path.join(dirout, "%s%s"%(fam, gftag)), 'w')
		fout.write('[FAMILIES]\n')
	# write out
	fout.write('- %s\n'%fam)
	fout.write('alignment = %s\n'%(nfali))
	fout.write('subst_model = %s\n'%(model))
	if nfgt: fout.write('starting_gene_tree = %s\n'%(nfgt))
	if nfmap: fout.write('mapping = %s\n'%(os.path.join(dirgt, bnmap)))
	if perfam:
		fout.close()

if not perfam:
	fout.close()
