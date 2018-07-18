#!/usr/bin/python
from parseALErec import getOrthologues
import ptg_utils as ptg
import tree2
import os, glob, sys

methods = ['strict', 'unicopy', 'mixed']
 
alerecdir = sys.argv[1]
outortdir = sys.argv[2]

for d in methods:
	pd = os.path.join(outortdir, d)
	if not os.path.isdir(pd):
		os.mkdir(pd)

foutdiffog = open(os.path.join(outortdir, 'diff_ortho_methods'), 'w')
foutdiffog.write('\t'.join(['family', 'nOG_strict', 'nOG_unicopy', 'nOG_mixed', 'overlap_strict_unicopy', 'overlap_strict_mixed', 'overlap_unicopy_mixed')+'\n')

#~ refspetree = tree2.AnnotatedNode(file='%s/ref_species_tree_41NeoPseudo_dated.nhx'%alerecdir, namesAsNum=True)

def summaryOGs(ogs, dlabs, N):
	n = len(ogs)
	print "number of OGs:", n
	cov = sum([len(x) for x in ogs])
	print "coverage of leaves:", cov, '/', N
	if cov != N:
		raise ValueError, "unclassified leaves:%s"%(repr(set(dlabs.values()) - set(reduce(lambda x,y: list(x)+list(y), ogs, []))))
	return n

lnfrec = glob.glob('%s/*ale.ml_rec'%(alerecdir))
for nfrec in lnfrec:
	fam = os.path.basename(nfrec).split('-', 1)[0]
	
	# just sample the first tree
	with open(nfrec, 'r') as frec:
		for i in range(12): frec.readline()
		recgtline = frec.readline()
		recgenetree = tree2.AnnotatedNode(nwk=recgtline, namesAsNum=True)
	#~ spetree, subspetree, lrecgt, recgtlines, restrictlabs, dnodeevt = parseALERecFile(nfrec, skipEventFreq=True)
		
	print "\n# # # #\n"
	N = recgenetree.nb_leaves()
	print "strict_ogs"
	strict_ogs, dlabs1 = getOrthologues(recgenetree, ALEmodel='dated', method='strict', verbose=True)
	n1 = summaryOGs(strict_ogs, dlabs1, N)
	print "\n# # # #\n"
	print "unicopy_ogs"
	unicopy_ogs, dlabs2 = getOrthologues(recgenetree, ALEmodel='dated', method='unicopy', verbose=True) #, refspetree=refspetree
	n2 = summaryOGs(unicopy_ogs, dlabs2, N)
	print "\n# # # #\n"
	print "mixed_ogs"
	mixed_ogs, dlabs3 = getOrthologues(recgenetree, ALEmodel='dated', method='mixed', verbose=True) #, refspetree=refspetree
	n3 = summaryOGs(mixed_ogs, dlabs3, N)
	o12 = sum([int(o in strict_ogs) for o in unicopy_ogs])
	o13 = sum([int(o in strict_ogs) for o in mixed_ogs])
	o23 = sum([int(o in unicopy_ogs) for o in mixed_ogs])
	print "overlap strict_ogs with unicopy_ogs:", o12
	print "overlap strict_ogs with mixed_ogs:", o13
	print "overlap unicopy_ogs with mixed_ogs:", o23
	
	foutdiffog.write('\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n'%(fam, n1, n2, n3, o12, o13, o23))
	
	for node in recgenetree:
		if node.is_leaf():
			node.edit_label(dlabs1[node.label()])
		else:
			node.edit_label('')
	dogs = {'strict':strict_ogs, 'unicopy':unicopy_ogs, 'mixed':mixed_ogs}
	for method, ogs in dogs.iteritems():
		ptg.colour_tree_with_leaf_groups(recgenetree, ogs)
		#~ ptg.colour_tree_with_constrained_clades(recgenetree, ogs)
		recgenetree.write_nexus(nftree.replace(".nhx", "_%s_orthologous_groups.nex"%method), ignoreBS=True)
		with open(os.path.join(outortdir, method, "%s_%s.orthologs"%(fam, method), 'w') as foutort:
			foutort.write('\n'.join([' '.join(x) for x in ogs])+'\n')
	
	print "\n# # # # # # # #\n # # # # # # # \n"

foutdiffog.close()
