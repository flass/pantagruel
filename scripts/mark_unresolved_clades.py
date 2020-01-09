#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""Select tree clades based on branch/node attributes and define constraint clades for MrBayesfor further bayesian tree inference, 
or produces alignments with collapsed clades (i.e. replace sets of sequences with a single representative of the clade).
"""

import os, sys, glob, getopt
import tree2, copy
import multiprocessing
from Bio import AlignIO, Align, Alphabet
#~ from random import shuffle
from random import randint
from ptg_utils import *

sys.setrecursionlimit(20000)

# a map of depth of current big recursive objects (trees) worked upon by each process
currentmaxreccursdepths = {} 
###### !!!! is this really sahred by mp processes ???? for the moment one when main(..., isparallel=False)

# functions to edit and monitor it, and to set Python maximum recursion limit accordingly
def checkmaxdepthandsetreclim(depths, mindefault=10000):
	maxval = max(depths.values()+[mindefault])
	if sys.getrecursionlimit() != maxval:
		sys.setrecursionlimit(maxval)
		print "# sys.setrecursionlimit(%d)"%maxval

def adddepth(depths, key, val):
	depths[key] = val
	checkmaxdepthandsetreclim(depths)

def rmdepth(depths, key):
	del depths[key]
	checkmaxdepthandsetreclim(depths)
# output folders and file extension for main function
doutdext = {'mbc':('mbconstraints', 'mbconstraints'), \
			'ccs':('constrained_clade_subalns', 'outgroup'), \
			'rle':('representative_leaves', 'representative_leaves'), \
			'col':('collapsed_alns', 'collapsed'), \
			'cgt':('coloured_genetrees', 'full_genetree.coloured'), \
			'nwk':('collapsed_ML_genetrees', 'collapsed')
		   }

##core functions
	
def select_clades_on_conditions(tree, clade_stem_conds, within_clade_conds, depth=1, testLeaves=True, testRoot=False, pruneSelected=True, nested=False, inclusive=True, order=2, verbose=2):
	"""select tree nodes based on generic conditions and return corresponding leaf sets.
	
	The conditions for selecting a clade are separated in two types relating to the the stem and within the clade, respectively.
	For each, a list of condition tuple(s) must be given, having the following structure: 
	
	clade_stem_conds=[(crit, operator, threshold_val [, leaf_mul]), ...]
	within_clade_conds=[(groupfun, crit, operator, threshold_val, depth), ...]
	
	where:
	- 'crit' is the name of a node/nbranch attribute, 
	- 'threshold_val' is the limit value, 
	- 'leaf_mul' is an OPTIONAL multiplicator to be applied to this threshold when the node is a leaf (default to 1.0)
	- 'operator', one of the following strings: '<', '>', '<=', '>=', '==',
	  indicating the logical test to do between the two former variables,
	- 'groupfun' is one of the following strings: 'mean', 'min', 'max', 
	  indicating the aggregating function for the several criterion values from within the clade,
	- 'depth' is the maximum rank of nodes below the clade stem from which 
	  the criterion values are gathered (negative depth means all nodes below).
	
	Note ALL conditions are to be met, i.e. they are combined with a logical 'AND'.
	
	Example:
		clade_stem_conds=[('bs', '>=', 80, 1)]
		within_clade_conds=[('max', 'lg', '<', 0.0001, -1), ('min', 'nb_leaves', '>=', 2, 1)]
	
	When 'nested' is True, no constraint selected nodes can be below another, 
	and resulting clades/leaf sets can be nested/included within the other node's.
	
	If 'pruneSelected' is True, the node that were selcted are removed from the observed tree.
	In a binary tree, this leads to fusion of the sister branch with the parent node's branch,
	which can impact the clade evaluation to follow - typically the sister clade will appear
	with a longer stem branch.
	
	
	If 'nested' is True and 'inclusive' is False, populations can be nested, but not overlap 
	- the nesting population will thus be paraphyletic; 
	if 'nested' is True and 'inclusive' is True, all populations remain monophyletic and can overlap.
	
	'order' defines he order of the tree traversal, which impacts how above rules ('nested' and 'inclusive') 
	are applied; the default order=2 means post-order traversal, thus leaves explored before their ancestors. 
	Other settings are not expected to yield anything desirable.
	"""
	if order!=2: print "Warning [mark_unresolved_clades.select_clades_on_conditions]: other tree traversal setting than default (2: post-order traversal) is not expected to yield anything desirable."
	
	def stem_testfun(node):
		"""generic test function for selection condition at the clade stem"""
		#~ print repr(clade_stem_conds)
		for clade_stem_cond in clade_stem_conds:
			try:
				crit, ope, thresh = clade_stem_cond
				leafmul = 1.0
			except ValueError:
				# case where optional leaf_mul is supplied as 4th value
				crit, ope, thresh, leafmul = clade_stem_cond
			if node.is_leaf():
				if crit in ['bs', 'boot']:
					# leaf branch support is always supperior than threshold (emulates maximal value)
					critval = thresh+1
				else:
					critval = getattr(node, crit)
					thresh = float(thresh) * float(leafmul)
			else:
				critval = getattr(node, crit)
			if callable(critval): val = critval()
			else: val = critval
			teststr = "%s %s %s"%(str(val), ope, str(thresh))
			etest = eval(teststr)
			if verbose>1: print "stem", node.label(), ':', teststr, etest
			if not etest: return False
		else:
			return True
			
	def within_testfun(node, constraintnodes):
		"""generic test function for selection condition whithin the clade"""
		if node.is_leaf():
			return True
		else:
			for grpfun, crit, ope, thresh, dp in within_clade_conds:
				lv = (None if (crit=='bs') else 'same')
				lsubval = node.getattr_down_n_nodes(crit, dp, ommitself=True, leafval=lv, stopatnodes=constraintnodes)
				if not lsubval: return False
				grpval = eval(grpfun)(lsubval)
				if grpval is None: return False
				llabs = node.get_leaf_labels()
				teststr = "%s %s %s"%(str(grpval), ope, str(thresh))
				etest = eval(teststr)
				if verbose>1: print "within", node.label(), ':', teststr, etest
				if not etest: return False
			else:
				return True
	
	def select_clades(tree):
		t = copy.deepcopy(tree)
		t.complete_internal_labels(order=0, ffel=True) # pre-order traversal numeration
		selectedclades = []
		selectedleaves = set([])
		selectednodes = []
		# iterate over all clades present in the tree (but the root)
		allnodes = t.get_sorted_children(order=order)
		n = 0
		#~ for node in allnodes:
		while n < len(allnodes):
			# from tips to root, accounting for the nodes that were removed
			node = allnodes[n]
			#~ print 'node', node
			n += 1 # for next iteration
			if (node.is_leaf() and not testLeaves) or (node.is_root() and not testRoot):
				continue # the while node loop
			if stem_testfun(node):
				# clade meets stem condition, e.g. it is enough supported to be constrained,
				# so can test the within-clade conditions (e.g. all-low branch supports)
				if nested: wt = within_testfun(node, selectednodes)
				else: wt = within_testfun(node, [])
				if wt:
					if not nested:
						for cnode in node.get_sorted_children(order=order):
							if cnode in selectednodes:
								# this clade would be nesting previously defines ones, skip
								continue # the for node loop
					leaves = node.get_leaf_labels()
					if verbose: print ' -> CC!', leaves
					if not inclusive:
						# restrict the leaf set to those not yet covered by selected leaf sets
						leaves = list(set(leaves) - selectedleaves)
						if not leaves:
							# all leaves are assigned to nested clades
							continue # the while node loop
					selectedclades.append(leaves)
					if pruneSelected: 
						# father will disapear from the node list
						if verbose: print 'pop', node.label(), node.newick(ignoreBS=1)
						collapsenodelabs = t.pop(node, tellReplacingNode=True)
						rmnode = t[collapsenodelabs[0]] 
						if verbose:
							print 'collapsenodelabs', collapsenodelabs, 
							print '; rmnode', rmnode.newick(ignoreBS=True)
						allnodes.remove(rmnode)
						if testRoot and (node is t):
							# cannot prune the root; rather return an empty node
							t = t.newnode()
						else:
							t.pop(node)
						# it takes time evaluating node nesting; no need to do it if pruning the node
						# so do not store selcted nodes/leaf sets
					else:
						selectedleaves |= set(leaves)
						selectednodes.append(node)
		return (t, selectedclades)
	
	if verbose: print "'select_clades_on_conditions(clade_stem_conds='%s', within_clade_conds='%s')"%(clade_stem_conds, within_clade_conds)
	tr = tree
	selectedclades = []
	nl0 = nl = tr.nb_leaves()
	nl1 = nl0+1
	nround = 0
	if verbose: print "Starts with %d leaves in the tree"%nl
	while nl1 > nl:
		nround += 1
		if verbose: print '# round %d ...'%nround
		nl1 = tr.nb_leaves()
		# may need two (or more?) passes to cover clades in the version with pruning
		tr, sc = select_clades(tr)
		nl = tr.nb_leaves()
		selectedclades += sc
		nspeinsc = sum([len(c) for c in selectedclades])
		if verbose:
			print '# round %d ... complete'%nround
			print 'leaves in selected clades:', nspeinsc
			print 'number of selected clades:', len(selectedclades)
			print 'remaining tree (%d leaves):\n'%nl, tr.newick(ignoreBS=True, comment='boot')
	else:
		nl = tr.nb_leaves()
		if nl==1 and not tr.get_leaf_labels()[0]:
			# empty tree
			nl = 0
	if pruneSelected:
		if (nspeinsc+nl)!=nl0:
			raise ValueError, "different counts of species in clades recorded: 'nspeinsc'+'nl': %d+%d; 'nl0': %d"%(nspeinsc, nl, nl0)
	return selectedclades
	
def mark_unresolved_clades(tree, cladesupport=None, subcladesupport=None, crit='bs', clade_stem_conds=None, within_clade_conds=None, nested=True, inclusive=False, pruneSelected=False, depth=1, verbose=2, **kw):
	"""Mark well supported clades with no strong support for the topology within the clade.
	
	i.e. the support for the clade is higher than 'cladesupport' and the supports 
	for nodes that are below this clade's node are less than the threshold 'minsubcladesupport'.
	
	This function is a wrapper for select_clades_on_conditions().
	
	The conditions for selecting he clade can be generalized to any criterion, 
	either by specifying 'crit' ('l' for branch lengths, 'boot' for branch support, ...)
	or by passing a list of condition tuples to 'clade_stem_conds' and 'within_clade_conds' arguments.
	(cf. select_clades_on_conditions() function)
	
	Working based on the consensus tree topology is justfied as long as 'cladesupport' >= 0.5,
	as all clades with PP > 0.5 will be represented in the consensus tree
	
	When 'nested' is False, no constraint clade/population can be nested or included in another.
	
	If 'nested' is True and 'inclusive' is False, populations can be nested, but not overlap 
	- the nesting population will thus be paraphyletic; 
	if 'nested' is True and 'inclusive' is True, all populations remain monophyletic and can overlap.
	
	'depth' define the maximum rank of nodes to explore below the clade's node 
	to evaluate low supports. For 'depth' >= 0, that many ranks will be explored; 
	if depth < 0, all splits below are explored, but those already alocated to other populations.
	'depth' defaults to 1, that is to just check resolved topology immediately below the node.
	('depth' argument is overridden by the value of fully-specified condition passed via the 'within_clade_conds' argument)
	"""
	# define generic test conditions for selection at the clade stem
	if clade_stem_conds is None:
		if cladesupport is None: raise ValueError, "must specify at least one of 'cladesupport' or 'clade_stem_conds' arguments"
		csc = [(crit, '>=', cladesupport)]
	else:
		if isinstance(clade_stem_conds, str): csc = eval(clade_stem_conds)
		else: csc = clade_stem_conds
	# define generic test conditions for selection whithin the clade
	if within_clade_conds is None:
		if subcladesupport is None: raise ValueError, "must specify at least one of 'subcladesupport' or 'within_clade_conds' arguments"
		wcc = [('max', crit, '<', subcladesupport, depth)]
	else:
		if isinstance(within_clade_conds, str): wcc = eval(within_clade_conds)
		else: wcc = within_clade_conds
	print 'csc:', csc, 'wcc:', wcc
	return select_clades_on_conditions(tree, clade_stem_conds=csc, within_clade_conds=wcc, nested=nested, inclusive=inclusive, pruneSelected=pruneSelected, testLeaves=False, verbose=verbose)

### output functions
		
def restrict_alignment_representative_leaves(constraints, tree, nffullali, dirout, radout, selectRepr=0, aliformatin='nexus', aliformatout='nexus', verbose=False, doutdext=doutdext, **kw):
	verbose = kw.get('verbose')
	supressout = kw.get('supressout', {})
	didseq = kw.get('didseq', {})
	aln = AlignIO.read(nffullali, format=aliformatin, alphabet=Alphabet.generic_dna)
	if verbose: print os.path.basename(nffullali)
	lcollapsedseqalnrowids = []
	loutgroups = []
	representativeseqalnrowids = []
	for i, constraint in enumerate(constraints):
		cladename = "clade%d"%i
		# collect clade sequences
		collapsedseqalnrowids = findSeqRecordIndexesFromSeqNames(aln, constraint)
		lcollapsedseqalnrowids.append(collapsedseqalnrowids)
		alicollapsedseqs = Align.MultipleSeqAlignment([seq for k,seq in enumerate(aln) if k in collapsedseqalnrowids])
		# select outgroup sequence: first leaf of sister clade
		outgroupleaf = findCladeOutgroup(constraint, tree, didseq)
		loutgroups.append(outgroupleaf)
		if not 'ccs' in supressout:
			outgroupseqalnrowid = findSeqRecordIndexesFromSeqNames(aln, outgroupleaf)
			alicollapsedseqs.append(aln[outgroupseqalnrowid])
			ccsoutd, ccsext = doutdext['ccs']
			AlignIO.write(alicollapsedseqs, os.path.join(dirout, ccsoutd, radout+"-%s-%s-%s.%s"%(cladename, ccsext, outgroupleaf, dalnext[aliformatout])), aliformatout)
		# select representative sequence out of the clade
		if callable(selectRepr): reprseq = selectRepr(alicollapsedseqs) # function returning a seq
		else: reprseq = alicollapsedseqs[selectRepr]
		reprseqalnrowid = findSeqRecordIndexesFromSeqNames(aln, reprseq.id)
		representativeseqalnrowids.append(reprseqalnrowid)
		if verbose: 
			for thing in ['i', 'cladename', 'leaves', 'reprseq.id', 'reprseqalnrowid']: #, 'outgroupleaf', 'outgroupseqalnrowid'
				print thing+':', eval(thing),
				print ''
	allcollapsedbutreprseqalnrowids = reduce(lambda i,j: set(i) | set(j), lcollapsedseqalnrowids, set([])) - set(representativeseqalnrowids)
	representativeseqids = [aln[reprseqalnrowid].id for reprseqalnrowid in representativeseqalnrowids]
	if not 'rle' in supressout:
		rleoutd, rleext = doutdext['rle']
		with open(os.path.join(dirout, rleoutd, radout+'-'+rleext), 'w') as foutrepr:
			for i, reprseqalnrowid in enumerate(representativeseqalnrowids):
				cladename = "clade%d"%i
				if verbose: print cladename, reprseqalnrowid, representativeseqids[i]
				foutrepr.write('\t'.join([cladename, representativeseqids[i]])+'\n')
				aln[reprseqalnrowid].id = cladename	# rename the sequence in the main alignment
	if not 'col' in supressout:
		coloutd, colext = doutdext['col']
		aliremainingseqs = Align.MultipleSeqAlignment([seq for k,seq in enumerate(aln) if k not in allcollapsedbutreprseqalnrowids])
		AlignIO.write(aliremainingseqs, os.path.join(dirout, coloutd, radout+"-%s.%s"%(colext, dalnext[aliformatout])), aliformatout)
		if len(aliremainingseqs) < 4:
			sys.stderr.write( "Warning: collapsed alignment '%s' has too few remaining sequences (%d out of %d) for subsequent tree building; however further sequence may be added in the clade replacement step\n"%(os.path.basename(nffullali), len(aliremainingseqs), len(aln)) )
	return (loutgroups, representativeseqids)

def collapse_tree_from_constraint_list(tree, constraints, representativeseqids, verbose=False):
	# first copy the tree to collapse so to leave untouched the input tree object
	if verbose: print "generate collapse tree from constraint list"
	coltree = copy.deepcopy(tree)
	assert len(constraints) == len(representativeseqids)
	i = -1
	for constraint, reprseqid in zip(constraints, representativeseqids):
		i += 1
		cladename = "clade%d"%i
		if verbose: print cladename
		for seqid in constraint:
			node = coltree.labelgetnode(seqid)
			if not node:
				# one of the identical sequences not present in the input tree
				if verbose: print "  %s not present in gene tree"%seqid
				continue # the for seqid loop
			if seqid == reprseqid:
				# make this leaf the representative of the clade
				if verbose:
					print "  %s is CC representative sequence; change label:"%seqid, node, '->', 
				node.edit_label(cladename)
				if verbose: print node
			else:
				if verbose: print "  %s in CC: pop leaf"%seqid
				coltree.pop(node)
	return coltree

def write_out_MrBayes_clade_constraints(constraints, outgroup, nfout, allowin=False, verbose=False, ilist=None, **kw):
	if isinstance(outgroup, list): outgroupleaves = outgroup
	elif isinstance(outgroup, str): outgroupleaves = [outgroup]
	else: raise ValueError, "unapropriate type for outgroup, must be a (list of) sequnce labels (str): %s"%repr(outgroup)
	# pick/verify outgroup leaf (is) out of constrained clades
	allleavesinconstraints = reduce(lambda i,j: set(i) | set(j), constraints, set([]))
	leftleaves = list(set(outgroupleaves) - allleavesinconstraints)
	if leftleaves: firstoutgroup = leftleaves[0]
	elif allowin: firstoutgroup = outgroupleaves[0]
	else: raise ValueError, "all outgroup sequence are located wihin a constraint clade"
	if verbose: print 'firstoutgroup:', firstoutgroup
	with open(nfout, 'w') as fout:
		if (firstoutgroup != 'None'): fout.write("outgroup %s\n"%(firstoutgroup))
		clades = []
		if isinstance(ilist, list): il = ilist
		else: il = range(len(constraints))
		for c, constraint in enumerate(constraints):
			i = il[c]
			fout.write("constraint clade%d -1 = %s\n"%(i, ' '.join(constraint)))
			clades.append("clade%d"%i)
		fout.write("prset topologypr=constraints(%s)\n"%(', '.join(clades)))

def main(nfgenetree, diraln, dirout, outtag, mkdircons=True, **kw):
	aliformatin = kw.get('aliformatin')
	diridentseq = kw.get('diridentseq')
	isparallel = kw.get('isparallel')
	print nfgenetree
	bnspl = os.path.basename(nfgenetree).split('.')
	if bnspl[0].startswith('RAxML_'): bngt = bnspl[1]
	else: bngt = bnspl[0]
	globaln = "%s/%s*%s*"%(diraln, bngt, aliformatin[:3])
	try:
		nfaln = glob.glob(globaln)[0]
	except IndexError:
		globaln = "%s/%s*aln*"%(diraln, bngt)
		nfaln = glob.glob(globaln)[0]
	print nfaln
	bnfaln = os.path.basename(nfaln).split('.')[0]
	if bnspl[0].startswith('RAxML_rootedTree'):
		# tree is already rooted, but the branch supports are storred in the comments
		genetree = tree2.AnnotatedNode(file=nfgenetree, keep_comments=True)
		for n in genetree:
			if str(n.comment()).isdigit():
				n.set_bs(float(n.comment()))
	else:
		 # trees is unrooted
		genetree = tree2.read_check_newick(nfgenetree, treeclass='AnnotatedNode')
		# tree is interpreted here as trifurcated at the root ; root it.
		genetree.resolveNode(outgroups='subroot')
	# there will be deepcopy operation on the tree, either to save its state before pruning (pop) below or in select_clades.
	# a deepcopy operation on a recursive tree2.Node object induces a cycle of ~7 function calls per nested node
	# knowing that there are (2*n)-1 nodes in a tree (n being the nuber of tree leaves),
	# one should set the recursion limit >> 7*2*n ; on sequential calls, set it to 10*2*n to be on the safe side with overheads of higher level function calls.
	if not isparallel: adddepth(currentmaxreccursdepths, bngt, 10*2*genetree.nb_leaves())
	# deal with potential information on sets of identical sequences 
	didseq = {}
	if diridentseq:
		# parse pairs of (reference, redundant) sequences that are identical
		globidseq = "%s/*%s*"%(diridentseq, bngt)
		gnfidseq = glob.glob(globidseq)
		if not gnfidseq: raise OSError, "cannot find file matching pattern: '%s'"%globidseq
		nfidseq = gnfidseq[0]
		with open(nfidseq, 'r') as fidseq:
			for line in fidseq:
				refidseq, redidseq = line.rstrip('\n').split('\t')
				didseq.setdefault(refidseq, []).append(redidseq)
		if didseq:
			# remove any redundant sequence from the gene tree before processing
			cleangenetree = copy.deepcopy(genetree)
			gtleaves = set(cleangenetree.get_leaf_labels())
			for refidseq, redidseqs in didseq.iteritems():
				for redidseq in redidseqs:
					if redidseq in gtleaves:
						cleangenetree.pop(redidseq)
						gtleaves.remove(redidseq)
		else:
			cleangenetree = genetree	
	else:
		cleangenetree = genetree	
	# detect unresolved clades
	constraintswithsingles = mark_unresolved_clades(cleangenetree, **kw) #, pruneSelected=True, inclusive=True
	if verbose: print 'constraintswithsingles =', constraintswithsingles
	# add to identical sequence map to the constrained clades definitions
	newconstraintsfromidseqs = []
	for refidseq, redidseqs in didseq.iteritems():
		# scan for existing clade that would contain the reference
		for c in constraintswithsingles:
			if refidseq in c:
				c += redidseqs
				break
		else:
			newconstraintsfromidseqs.append([refidseq]+redidseqs)
	if verbose: print 'newconstraintsfromidseqs =', newconstraintsfromidseqs
	# for reporting, filter out contraint clades that are just made of one leaf (NB: these are useful for proper definition of other constraint clades, when nested, non-inclusive clades are allowed)
	constraints = [c for c in constraintswithsingles+newconstraintsfromidseqs if len(c)>1]
	if verbose: print 'constraints =', constraints
	# write out subalignments and the main alignment with collapsed clades
	loutgroups, representativeseqids = restrict_alignment_representative_leaves(constraints, genetree, nfaln, dirout, radout=bnfaln, selectRepr=0, didseq=didseq, **kw)
	if not 'mbc' in supressout:	
		mbcoutd, mbcext = doutdext['mbc']
		if mkdircons is True:
			mbcoutd = os.path.join(mbcoutd, bnfaln)
			dout = os.path.join(dirout, mbcoutd)
			if not os.path.isdir(dout):
				os.mkdir(dout)
		for i, constraint in enumerate(constraints):
			cladename = "clade%d"%i
			# write out MrBayes clade constraint for the sub-alignment, in order to compute subalignment samples and/or ancestral sequence
			write_out_MrBayes_clade_constraints([constraint], loutgroups[i], os.path.join(dirout, mbcoutd, bnfaln+'-'+cladename+'.'+mbcext), ilist=[i], verbose=verbose)
	if not 'nwk' in supressout:
		colgenetree = collapse_tree_from_constraint_list(genetree, constraints, representativeseqids, verbose=verbose)
		nwkoutd, nwkext = doutdext['nwk']
		colgenetree.write_newick(os.path.join(dirout, nwkoutd, bnfaln+'-%s.nwk'%nwkext), ignoreBS=True)
	if not 'cgt' in supressout:	
		fmtcoltree = kw.get('format_color_tree')
		cgtoutd, cgtext = doutdext['cgt']
		colour_tree_with_constrained_clades(genetree, constraints, force=True)
		genetree.complete_internal_labels(order=0, ffel=True)
		if fmtcoltree.lower() in ['xml', 'phyloxml']:
			genetree.write_phyloXML(os.path.join(dirout, cgtoutd, bnfaln+'-%s.xml'%cgtext), ignoreBS=True)
		elif fmtcoltree.lower() in ['nex', 'nexus']:
			genetree.write_nexus(os.path.join(dirout, cgtoutd, bnfaln+'-%s.nex'%cgtext), ignoreBS=True)
		else:
			raise ValueError, "specified format '%s' for output coloured-branch tree is not valid; please select among '[phylo]xml' or 'nex[us]'"%fmtcoltree
	# done risking going over reccursion limt
	if not isparallel: rmdepth(currentmaxreccursdepths, bngt)
		
	
def usage():
	s = 'Usage:\n'
	s += '\tpython mark_unresolved_clades.py --in_gene_tree_list=/path/to/genetree_filepath_list [OPTIONS]\n'
	s += '\n  Options:\t\t[Value:]\n'
	s += '\n    1. general options:\n'
	s += '\n      a. input:\n\n'
	s += '  --in_gene_tree_list\tpath\tpath to file containing the list of paths for source gene trees (MANDATORY).\n'
	s += '\t\t\t\t (which must be point estimates [ML, consensus] in Newick format).\n'
	s += '  --diraln\t\tpath\tpath to source alignment folder; defaults to the same older as each listed gene tree.\n'
	s += '  --dir_identseq\tpath\tpath to folder containing lists of identical sequences in the gene family;\n'
	s += '\t\t\t\teach file must be formated as a tab delimited table with two columns, indicating reference and redundant sequences, respectively.\n'
	s += '\t\t\t\tIf the redundant sequences are present in the input gene tree, they will be removed before processsing.\n'
	s += '\t\t\t\tIn any case, the reference:(redundant1, redundant2, ...) sequence name map is later integrated with the map of clades to collapse.\n'
	s += '  --fmt_aln_in\t\tfmt\tformat of input alignments; defaults to \'nexus\'.\n'
	s += '\n      b. output:\n\n'
	s += '  --dirout\t\tpath\tspecify output folder; defaults to the same older as each listed gene tree.\n'
	s += '  --out_mbcons_tag\tstr\tsuffix for MrBayes constraint output files; defaults to \'.mbconstraints_%.2g_%.2g\'%(cladesupport, subcladesupport).\n'
	s += '  --fmt_aln_out\t\tfmt\tformat of output alignments; defaults to \'nexus\'\n'
	s += '  --fmt_col_genetrees\tfmt\tspecify format {phyloxml|nexus} for a version of\n'
	s += '\t\t\t\tthe (full) input gene tree with constrained clades coloured; default to \'nexus\'.\n'
	s += '  --no_%x_output\t\tsuppress writing the corresponding \'dirout/%x\' output folder, with %x in:\n'
	s += '\t\t\t\t {%s}.\n'%(', '.join([outd for outd, ext in doutdext.values()]))
	s += '\n      c. misc:\n\n'
	s += '  --threads\t\tint\tnumber of parallel threads; defaults to 1\n'
	s += '  --verbose\t\t\tset verbose mode on.\n'
	s += '  --help\t\t\tprint this mesage and exit.\n'
	s += '\n    2. to select unresolved clades:\n'
	s += '\n      a. simple way; combine both \'--clade_support\' and \'--subclade_support\' options (and falcultatively \'--criterion\' and \'--depth\'):\n\n'
	s += '  --clade_support\tfloat\tdefine the minimum clade support (on the clade\'s stem branch); default with value of 0.5\n'
	s += '  --subclade_support\tfloat\tdefine the maximum branch support under the clade\'s stem (over d branches, see \'--depth\'); default with value of 0.2\n'
	s += '  --criterion\t\tstr\tdefine which criteion is evaluated, default to \'bs\' for branch support; can be \'lg\' for branch length, etc.\n'
	s += '  --depth\t\tint\tdefine depth within the clades subtree over which the criterion is evaluated\n'
	s += '\n      b. more formal way; describe completely the conditions at the stem and whitin the clade; allows to combine different criteria:\n\n'
	s += '  --clade_stem_conds\tstr\tconditions on clade stem (python-interpretable str); overrides values from \'--clade_support\'.\n'
	s += '  --within_clade_conds\tstr\tconditions on clade subtree (python-interpretable str); overrides values from \'--subclade_support\' and \'--depth\'.\n'
	s += '    Example:\n'
	s += '    	clade_stem_conds="[(\'bs\', \'>=\', 80)]"\n'
	s += '    	within_clade_conds="[(\'max\', \'lg\', \'<\', 0.0001, -1), (\'min\', \'nb_leaves\', \'>=\', 2, 1)]"\n'
	return s
	
if __name__=='__main__':

	supressoutopt = ['_'.join(['no', outd, 'output']) for outd, ext in doutdext.values()]
	
	opts, args = getopt.getopt(sys.argv[1:], 'hv:', ['in_gene_tree_list=', 'dir_identseq=', \
													 'clade_support=', 'subclade_support=', 'criterion=', 'depth=', \
													 'clade_stem_conds=', 'within_clade_conds=', \
	                                                 'diraln=', 'fmt_aln_in=', 'fmt_aln_out=', 'fmt_col_genetrees=', \
	                                                 'dirout=', 'out_mbcons_tag=', \
	                                                 'threads=', 'verbose=', 'help']+supressoutopt)
	dopt = dict(opts)
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	nflnfgenetree = dopt['--in_gene_tree_list']
	cladesupport = float(dopt.get('--clade_support', 0.5))
	subcladesupport = float(dopt.get('--subclade_support', 0.2))
	crit = dopt.get('--criterion', 'bs')
	depth = int(dopt.get('--criterion', 1))
	clade_stem_conds = dopt.get('--clade_stem_conds')
	within_clade_conds = dopt.get('--within_clade_conds')
	nthreads = int(dopt.get('--threads', 1))
	verbose = dopt.get('-v', dopt.get('--verbose', 0))
	outtag = dopt.get('--out_mbcons_tag', '.%s_stem%.2g_within%.2g'%(crit, cladesupport, subcladesupport))
	dirout = dopt.get('--dirout') 
	diraln = dopt.get('--diraln')
	alnfmtin = dopt.get('--fmt_aln_in', 'nexus')
	alnfmtout = dopt.get('--fmt_aln_out', 'nexus')
	fmtcoltree = dopt.get('--fmt_col_genetrees', 'nex')
	diridentseq = dopt.get('--dir_identseq')
	
	if clade_stem_conds and verbose:
		print "clade_stem_conds = %s"%repr(clade_stem_conds)
	if within_clade_conds and verbose:
		print "within_clade_conds = %s"%repr(within_clade_conds)
	
	# check accessibility of ouput dir
	for outd, ext in doutdext.values():
		if not os.path.isdir(os.path.join(dirout, outd)):
			os.mkdir(os.path.join(dirout, outd))
	
	supressout = []
	for key, val in doutdext.items():
		outd, ext = val
		if '_'.join(['--no', outd, 'output']) in dopt:
			supressout.append(key)
	
	# set the arguments to pass to pool.map() (parallel iteration)
	def mainSetArgs(nfgenetree):
		dgt = os.path.dirname(nfgenetree)
		dout = dirout if dirout else dgt
		daln = diraln if diraln else dgt
		main(nfgenetree, daln, dout, outtag, diridentseq=diridentseq, \
		     cladesupport=cladesupport, subcladesupport=subcladesupport, crit=crit, depth=depth, \
		     clade_stem_conds=clade_stem_conds, within_clade_conds=within_clade_conds, \
		     aliformatin=alnfmtin, aliformatout=alnfmtout, format_color_tree=fmtcoltree, \
		     supressout=supressout, verbose=verbose, isparallel=(nthreads > 1))
	
	with open(nflnfgenetree, 'r') as flnfgenetrees:
		lnfgenetrees = [line.strip(' \n') for line in flnfgenetrees]
		
	if nthreads > 1:
		pool = multiprocessing.Pool(processes=nthreads)
		pool.map(mainSetArgs, lnfgenetrees)
	else:
		for nfgenetree in lnfgenetrees:
			mainSetArgs(nfgenetree)
	
		
