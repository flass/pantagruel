#!/usr/bin/python

import tree2
import sys, getopt

def nodeAge(node, fromRoot=False):
	if fromRoot:
		return node.distance_root()
	else:
		oneleaf = node.get_one_leaf()
		return node.distance(oneleaf)

def filterNode(anchornode, excluding=[], excluding_clade=[]):
	lnodes = []
	for node in anchornode.get_children():
		nodeage = nodeAge(node, fromRoot=ageFromRoot)
		skipclade = (node in excluding_clade)
		skipnode  = ( (node in excluding) or skipclade \
		           or ((not (maxAge is None)) and (nodeage > maxAge)) \
		           or ((not (minAge is None)) and (nodeage < minAge)) )
		if not skipnode: lnodes.append(node)
		if not skipclade: lnodes += filterNode(node, excluding=excluding, excluding_clade=excluding_clade)
	return lnodes

opts, args = getopt.getopt(sys.argv[1:], 'hv', ['intree=', 'out=', \
                                                'older_than=', 'younger_than=', \
                                                'only_below=', 'excluding=',  'excluding_clade=', \
                                                'age_from_root', 'root_age=']

nfrefspetree = dopt['--intree']
nfout = dopt['--out']
minAge = dopt.get('--older_than')
maxAge = dopt.get('--younger_than')
onlyAncs = dopt.get('--only_below').split(',')
exclNodes = dopt.get('--excluding').split(',')
exclClades = dopt.get('--excluding_clade').split(',')
ageFromRoot = ('--age_from_root' in dopt) # if not (default), assumes that counting age from the leaves is valid, i.e. that the tree is ultrametric
rootAge = float(opt.get('--root_age', 0.0)) # if the above is set to true, this option defines the 

refspetree = tree2.AnnotatedNode(file=nfrefspetree)
refspetree.complete_internal_labels(order=0, ffel=True)
if not os.path.exists(os.path.dirname(nfout)):
	raise ValueError, "Directory for output file does not exists"

anchornodes = []
if onlyAncs:
	for onlyAnclab in onlyAncs:
		anchornodes.append(refspetree[onlyAnclab])
else:
	anchornodes.append(refspetree)

lnodes = []
for anchornode in anchornodes:
		lnodes += filterNode(anchornode, excluding=excluding, excluding_clade=excluding_clade)
		
with open(nfout, 'w') as fout:
	for node in lnodes:
		fout.write(node.label()+'\n')

refspetree.write_newick(nfout+'.reftree', ignoreBS=True)
