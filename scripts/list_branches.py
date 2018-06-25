#!/usr/bin/python

import tree2
import sys, getopt

def nodeAge(node, fromRoot=False):
	if fromRoot:
		return node.distance_root()
	else:
		oneleaf = node.get_one_leaf()
		return node.distance(oneleaf)

def filterNode(node, minAge=None, maxAge=None, excluding=[], excluding_clade=[], ageFromRoot=False):
	lnodes = []
	nodeage = nodeAge(node, fromRoot=ageFromRoot)
	skipclade = (node in excluding_clade)
	skipnode  = ( (node in excluding) or skipclade \
			   or ((not (maxAge is None)) and (nodeage > maxAge)) \
			   or ((not (minAge is None)) and (nodeage < minAge)) )
	if not skipnode: lnodes.append(node)
	if not skipclade:
		for cnode in node.get_children():
			lnodes += filterNode(cnode, minAge=minAge, maxAge=maxAge, excluding=excluding, excluding_clade=excluding_clade, ageFromRoot=ageFromRoot)
	return lnodes

opts, args = getopt.getopt(sys.argv[1:], 'hv', ['intree=', 'out=', \
                                                'older_than=', 'younger_than=', 'min_age=', 'max_age=', \
                                                'only_below=', 'excluding=',  'excluding_clade=', \
                                                'age_from_longest_tip', 'root_age='])

nfrefspetree = dopt['--intree']
nfout = dopt['--out']
minAge = dopt.get('--older_than', dopt.get('--min_age'))
maxAge = dopt.get('--younger_than', dopt.get('--max_age'))
onlyAncs = dopt.get('--only_below').split(',')
exclNodes = dopt.get('--excluding').split(',')
exclClades = dopt.get('--excluding_clade').split(',')
rootAge = float(opt.get('--root_age'))
ageFromMaxTip = ('--age_from_longest_tip' in dopt) # if not (default), assumes that counting age from any of the tips is equivalent,
                                                   # i.e. that the tree is ultrametric

if not os.path.exists(os.path.dirname(nfout)):
	raise ValueError, "Directory for output file does not exists"
refspetree = tree2.AnnotatedNode(file=nfrefspetree)
refspetree.complete_internal_labels(order=0, ffel=True)
if ageFromMaxTip:
	treeage = refspetree.max_leaf_distance()
else:
	treeage = refspetree.get_one_leaf().distance_root()
	print "Warning: will consider the given tree is ultrametric"
if not (rootAge is None):
	refspetree /= (treeage * rootAge) # scale the tree
else:
	print "Warning: the age of the tree is not scaled: root age is currently %f"%treeage

anchornodes = []
if onlyAncs:
	for onlyAnclab in onlyAncs:
		anchornodes.append(refspetree[onlyAnclab])
else:
	anchornodes.append(refspetree)

lnodes = []
for anchornode in anchornodes:
		lnodes += filterNode(anchornode, minAge=minAge, maxAge=maxAge, excluding=excluding, \
		                     excluding_clade=excluding_clade, ageFromRoot=ageFromMaxTip)

with open(nfout, 'w') as fout:
	for node in lnodes:
		fout.write(node.label()+'\n')

refspetree.write_newick(nfout+'.reftree', ignoreBS=True)
