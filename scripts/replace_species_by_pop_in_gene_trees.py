#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import tree2
import sys, os, getopt, glob
from ptg_utils import *
from mark_unresolved_clades import select_clades_on_conditions
from scipy import stats
import copy
import igraph

import multiprocessing as mp
from Bio import AlignIO
from Bio.Phylo import BaseTree
import traceback

phyloproftag='phyloprofiles'
maxreclimupfactor = 1 # range for tentative increase of the recursion limit that was set by default or via command-line argument when it is not enogh to read a gene tree; e.g. 1.5 allows for 50% increase, 1 allows for none

default_psc = [('lg', '>=', 0.0005), ('bs', '>=', 80)]
default_wpc = [('max', 'lg', '<=', 0.0005, -1)]

def speciesTreePopulations(spetree, pop_stem_conds, within_pop_conds, nested=False, inclusive=False, taglen=3, poptag=None, **kw):
	"""find supported clades in the species tree
	
	when inclusive is False, populations can be defined so as to emerges from another, with the more ancestral not be inclusive on of another, i.e. a population
	"""
	verbose = kw.get('verbose')
	lpops = select_clades_on_conditions(spetree , clade_stem_conds=pop_stem_conds, within_clade_conds=within_pop_conds, testRoot=True, nested=nested, inclusive=inclusive, **kw)
	# collect lone species not included in a population
	allspeinpops = set([])
	for pop in lpops: allspeinpops |= set(pop)
	for spe in spetree.get_leaf_labels():
		if not (spe in allspeinpops):
			if verbose: print "add lone species %s as a population"%spe
			lpops.append([spe])
	# generate pop names based on the major 3 first letter of species in it
	lnamepops = []
	dnamecount = {}
	for pop in lpops:
		assert (pop[0] is not None)
		if len(pop)>1:
			if poptag:
				maxtag = poptag
			else:
				taglet = [spe[:taglen] for spe in pop]
				tagcount = [(tag, taglet.count(tag)) for tag in set(taglet)]
				maxtag = sorted(tagcount, key=lambda x: x[1])[-1][0]
			dnamecount.setdefault(maxtag, 0)
			dnamecount[maxtag] += 1
			ntag = dnamecount[maxtag]
			# keep the opulation order, as potentially nested clades must be selected in the same order
			lnamepops.append(("%s-POP%d"%(maxtag, ntag), pop))
		else:
			lnamepops.append((pop[0], pop))
	return lnamepops

def collapsePopulationInSpeciesTree(spetree, lnamepops, fast=True, speciestoprune=[], nested=True, keepUltrametric=True, collapseAllPops=True, verbose=0):
	"""if population are nested, the order of the input population (name, (members,)) list matters!
	verbose=0: silent, verbose=1: general progression, verbose=2: detail of tree edition.
	"""
	if verbose:
		print 'collapsePopulationInSpeciesTree()'
		print "initiate with a species tree with %d leaves"%( spetree.nb_leaves())
	
	# first mark the population ancestor node
	spetree.complete_node_ids(order=2, force=True) # post-order traversal numeration
	poptree = copy.deepcopy(spetree) 
	#~ poptree = spetree.deepcopybelow(keep_lg=True)
	dnamepop = dict(lnamepops)
	# list of nodes independent of tree structure (i.e. pruning or unlinking nodes in the tree won't affect list sequence]
	lnamepopancs = [(popname, poptree.mrca(popspe)) for popname, popspe in lnamepops]
	if verbose: print 'lnamepopancs:\n%s'%( repr_long_list(lnamepopancs) )
		
	if speciestoprune:
		if verbose: print "start by pruning species from the tree:", speciestoprune #repr_long_list(list(speciestoprune))
		#~ try:
			#~ print speciestoprune
			#~ speciestoprune.sort(key=lambda x: poptree[x].nodeid())
		#~ except Exception, e:
			#~ print 'Caught exception (in thread %d):'%os.getpid()
			#~ traceback.print_exc()
			#~ sys.stdout.flush()
			#~ raise e
			
		for spetoprune in speciestoprune:
			s2p = poptree[spetoprune]
			if s2p:
				if verbose: print "prune %s: poptree.pop(%s)"%(spetoprune, repr(s2p)),
				popedspe = poptree.pop(s2p, verbose=max(verbose-1, 0))
				if verbose: print "pruned", popedspe
			else:
				raise IndexError, "species to prune '%s' is missing in poptree"%spetoprune
		else:
			if verbose: print ""
	
	if fast:
		# ensure sorting pop. ancestor nodes following post-order traversal numeration in the original tree
		lnamepopancs.sort(key=lambda x: x[1].nodeid())
		if verbose: print 'lnamepopancs (sorted according to popancs nodeids):\n%s'%( repr_long_list(lnamepopancs) )
	lpopname, lpopancs = zip(*lnamepopancs)
	lpopancids = [popanc.nodeid() for popanc in lpopancs] # increasing sequence of node ids
	
	for a, (popname, popanc) in enumerate(lnamepopancs):
		lpopspe = dnamepop[popname]
		p = popanc.nodeid()
		if verbose:
			print '#%d. %s (id: %d):'%(a, popname, p), lpopspe
			print ' in subtree:', popanc.newick(ignoreBS=True)
		if len(lpopspe)==1: 
			# leaf/single-species population; skip
			continue
		sublen = popanc.max_leaf_distance()
		if nested:
			if fast:
				# explore the presence of other pop. ancestor only within the subtree defined by the current node
				# which in post-order traversal is the slice of the node index list from the current node's first leaf up to itself (p-th position).
				# All nodes below a node (i.e. the node's subtree) being always before the node, and making an uninterupted series with only those nodes,
				# it follws that the slice is [(p-n+1):p], with n the number of nodes in the subtree (including the node itself, as returned by Node.nb_all_children()).
				# just have to translate that into the indexes of the pop. ancestor node id list,  which is a subset of the full node id set
				# <=> [k:a] with a = lpopancids.index(p) and k = min([t for t, j in enumerate(lpopancids) if j>=(p-n+1)])
				n = spetree.idgetnode(popanc.nodeid()).nb_all_children() # rely on parent-child structure of the original, uncut tree
				k = 0
				for k, j in enumerate(lpopancids):
					if j>=(p-n+1): break
				if verbose: print 'check subancs within lpopancs[%d:%d] i.e. from nodeid %d'%(k, a, j)
				subancs = popanc.is_parent_of_any(lpopancs[k:a], returnList=True)
			else:
				if verbose: print 'check all subancs'
				subancs = popanc.is_parent_of_any(lpopancs, returnList=True)
			if verbose: print 'subancs (labels):', [su.label() for su in subancs]
			if subancs:
				# remove population clade (or paraphyletic group) leaf by leaf
				# until only other population ancestor and one arbitrary leaf are left
				nkept = 0
				for j, spe in enumerate(lpopspe):
					if spe in popanc.children_labels():
						# leaf is just below the population ancestor 
						# must not pop() that leaf to avoid destroying the population ancestor node
						if nkept==0:
							# finding a leaf just under the population ancestor must happen at least once, for which the leaf is kept
							nkept += 1
							# mark the population name on that leaf; allow for extra branch length to be kept relative to the true population ancestor 
							# when non-nested populations have not, but that should be a small bias given condition defining populations should be based on 
							# the clade being a rake, i.e. all branch are of approx. null lentgh within the clade.
							# this also avoids creating a branch of actually zero length, which might lead to later trouble e.g. in ALE time slice definition.
							popanc[spe].edit_label(popname)
							if verbose: print 'kept', spe, 'as', popname
						else:
							# there is a multifurcation below the pop ancestor; simply unlink that child
							if verbose: print 'unlink', spe,
							popanc.unlink_child(popanc[spe], silent=1-(verbose-1))
						
					else:
						if j==(len(lpopspe)-1) and nkept==0:
							# finding a leaf just under the population ancestor has not yet happen,
							# but at least one leaf must be kept: mark the population name on that leaf.
							popanc[spe].edit_label(popname)
							if verbose: print 'kept', spe, 'as', popname
						else:
							# prune the leaf (using tree2.Node.pop()):
							# lead to deconnecting its father node from the rest of tree 
							# and connecting the brother node to the grand-father node
							if verbose: print 'pop', spe #,
							pspe = popanc.pop(spe)
							#~ print pspe
	
			else:
				# no other population ancestor below
				if set(popanc.get_leaf_labels())!=set(lpopspe):
					print 'set(popanc.get_leaf_labels()):'
					print set(popanc.get_leaf_labels())
					print 'set(lpopspe):'
					print set(lpopspe)
					print 'set(popanc.get_leaf_labels()) - set(lpopspe)'
					print set(popanc.get_leaf_labels()) - set(lpopspe)
					raise ValueError, "There should be no other population ancestor below '%s' pop ancestor node (in thread %d)"%(popname, os.getpid())
				popanc.as_leaf(newlabel=popname, silent=max(0, 1-(verbose-1)))
				if keepUltrametric and sublen!='NA': popanc += sublen	# reports distances to longest leaf lineage below the ancestor so that an ultrametric tree keeps being ultrametric
				if verbose: print 'as_leaf', popname, '+=', sublen
						
		else:
			# assume populations are monophyletic and the clade can be removed at once
			assert set(popanc.get_leaf_labels())==set(lpopspe)
			popanc.as_leaf(newlabel=popname, silent=1-(verbose-1))
	
	#~ print poptree
	if collapseAllPops and (not set(lpopname) == set(poptree.get_leaf_labels())):
		print "set(lpopname) - set(poptree.get_leaf_labels()):\n", set(lpopname) - set(poptree.get_leaf_labels())
		print "set(poptree.get_leaf_labels()) - set(lpopname):\n", set(poptree.get_leaf_labels()) - set(lpopname)
		raise IndexError, "population list and population tree are discordant"
			
	return poptree

def _popCountDict(lspe, dspe2pop, dnamepops, asList=True, addpopnames=[]):
	dspecount = {spe:lspe.count(spe) for spe in set(lspe)}
	dpopcount = {}
	if asList:
		# will return occurence profile (count per strain)
		popnames = set([dspe2pop[spe] for spe in dspecount]+addpopnames)
		for popname in popnames:
			dpopcount[popname] = [dspecount.get(spe, 0) for spe in dnamepops[popname]]
	else:
		# will return total occurence count
		for spe in dspecount:
			popname = dspe2pop[spe]
			dpopcount.setdefault(popname, 0) + dspecount[spe]
	return dpopcount

def combineMonophyleticPopulations(dpopcount, poptree, dnamepops, dpopnd={}, maxpopnodedist=3, **kw):
	"""try and pool represented populations if they are monophyletic (or paraphyletic but closely related) 
	and each have at least a moderately high occurence frequency of the gene (median >= 0.5).
	
	maxpopnodedist gives an extended criterion to consider the monophyly of populations, 
	based on node distance of clade ancestors:
	if maxpopnodedist < 2, NO populations will be lumped together;
	if maxpopnodedist = 2, the criterion is that only sister (i.e. truly monophyletic) populations will be lumped together;
	if maxpopnodedist > 2, it is possible that populations that are close in the tree but not sister can be joined;
	In any case, if of two populations, one ancestor is a descendant of the other ancestor, populations will be aggreagated.
	This is to account for imperfect gene conservation in clades, e.g. if in a clade of populations (((popA,popB),popC),popD),
	a gene is present (conserved) in popB, popC, and popD but absent in popA, suggesting the gene was ancestrally present 
	and lost in popA.
	"""
	verbose = kw.get('verbose', False)
	if maxpopnodedist < 2:
		if verbose: print "will not attempt to aggregate populations into monophyletic groups"
		return dpopcount 
	elif maxpopnodedist == 2:
		grouptype = 'monophyletic'
	else:
		grouptype = 'closely related'
	if len(dpopcount) > 1:
		pnames = [pname for pname in dpopcount if median(dpopcount[pname]) >= 0.5]
		# use a graph approach to cluster nodes in the tree that are under a threshold inter-node distance
		p2p = igraph.Graph()
		p2p.add_vertices(len(pnames))
		p2p.vs['name'] = pnames
		edges = []
		for i, popnamei in enumerate(pnames):
			ancpopi = poptree[popnamei]
			for j, popnamej in enumerate(pnames[i+1:]):
				pnd = dpopnd.get(popnamei, {}).get(popnamej)
				if pnd is None:
					ancpopj = poptree[popnamej]
					pnd = ancpopi.node_distance(ancpopj)
				if pnd <= maxpopnodedist:
					edges.append((popnamei, popnamej))
		p2p.add_edges(edges)
		compsp2p = p2p.components()
		for comp in compsp2p:
			if len(comp) > 1:
				if verbose: print "combine %s populations %s, repesented by a total of %f sequences"%(grouptype, repr(jointpopnames), float(jointpopcount))
				cpnames = tuple(p2p.vs[p]['name'] for p in comp)
				jointpopcount = reduce(lambda x,y: x+y, (dpopcount[cpname] for cpname in cpnames), [])
				dpopcount[cpnames] = jointpopcount
				for cpname in cpnames:
					del dpopcount[cpname]		
		#~ for i, popnamei in enumerate(pnames):
			#~ indcountsi = dpopcount[popnamei]
			#~ if median(indcountsi) >= 0.5:
				#~ jointpopcount = indcountsi
				#~ ancpopi = poptree[popnamei]
				#~ jointpopnames = addStrtoTuple(popnamei, ())
				#~ for j, popnamej in enumerate(pnames[i+1:]):
					#~ indcountsj = dpopcount[popnamej]
					#~ if median(indcountsj) >= 0.5:
						#~ # at least hlaf of the joint sample has per-strain occurence greater than 0
						#~ ancpopj = poptree[popnamej]
						#~ jpn = addStrtoTuple(popnamej, jointpopnames)
						#~ # based on full species tree
						#~ # jointleaves = reduce(lambda x,y: x+y, (dnamepops[pn] for pn in jpn))
						#~ # if ancpopj.is_child(ancpopi) or ancpopi.is_child(ancpopj): mono = True
						#~ # elif maxpopnodedist>2: mono = (dpopnd.get(popnamei, {}).get(popnamej, ancpopi.node_distance(ancpopj)) <= maxpopnodedist)
						#~ # elif spetree.is_monophyletic(jointleaves): mono = True
						#~ # more efficiently based on collapsed species tree = population tree
						#~ if maxpopnodedist>2: mono = (dpopnd.get(popnamei, {}).get(popnamej, ancpopi.node_distance(ancpopj)) <= maxpopnodedist)
						#~ else: mono = False
						#~ if mono:
							#~ jointpopnames = jpn
							#~ del dpopcount[popnamei]
							#~ del dpopcount[popnamej]
							#~ jointpopcount += indcountsj
							#~ dpopcount[jointpopnames] = jointpopcount
							#~ if verbose: print "combine populations %s"%repr(jointpopnames)
							#~ # reccursive call allow to iteratively join monophyletic group
							#~ dpopcount = combineMonophyleticPopulations(dpopcount, poptree, dnamepops, dpopnd=dpopnd, maxpopnodedist=maxpopnodedist, **kw)
							#~ return dpopcount
	return dpopcount

def getCladeAncestor(dpopcount, dnamepops, spetree, dspe2pop, medpopcountthresh=0.5, **kw):
	"""rank populations or population groups based on a sequence of criteria tofind the most likely population of origin of the sequences in the collapsed clade.
	
	determine what population (or group of populations) was the original donor of the sequence based on:
	# 1) the median gene occurence frequency being >= 0.5 AND the population census being > 1 (i.e. lone strains don't count) [True rank first] 
	# 2) the higher depth in the species tree of most recent common ancestor of the population(s) (= minimum cumulated branch distance from the root) [best with ultrametric] 
	# 3) the larger census (= count of species classified in the population(s), regarless of the gene presence/absence = sample size when estimating freq)
	# 4) the mean gene occurence frequency
	
	High frequncy in lone strain is disregarded in first criterion 
	because lone strain lineage tend to branch higher in the species tree, 
	and thus to always rank first on the 2nd criterion.
	riterion based on mean frequency is considered last only to break ties,
	because mean is too sensitive to sample size
	"""
	def eval_criteria(x, medth=medpopcountthresh):
		"""expect vector with (pop_name, popcount, pop_anc_node_in_spe_tree)"""
		medc = median(x[1]) ; meanc = mean(x[1]) ; lenc = len(x[1])
		return ((medc>=medth and lenc>1), -1*(x[2].distance_root()), lenc, meanc)
	
	verbose = kw.get('verbose', False)
	allpopdistr = [(popname, popcount, spetree[popname]) for popname, popcount in dpopcount.iteritems()]
	if len(allpopdistr) > 1:			
		## rank populations according to 
		# from best values (index 0) to worst (index -1), breaking ties with sample size (more = better) and mean freq (higher = better)
		allpopdistr.sort(key=lambda x: eval_criteria(x), reverse=True) 
		if verbose:
			print "\tancestor_to_root_dist\toccurence/sample_size\tpopname"
			for x in allpopdistr[:5]:
				print "\t".join(str(y) for y in ['', x[2].distance_root(), "%d/%d"%(sum(x[1]),len(x[1])), repr(x[0])])
			u = len(allpopdistr)-5
			if u>0: print "\t... (%d other pops)"%u
	ancpopname, indcounts, ancpopclade = allpopdistr[0]
	# from a proposed set of ancestral population, record the related populations that are missing the gene
	# but are located below the the MRCA of ancestral populations (putative losses since the MRCA)
	#~ lostpopnames = list((set(ancpopclade.get_children_labels()) & set(dnamepops.keys())) - set(addStrtoTuple(ancpopname, ())))
	lostpopnames = list(set(dspe2pop[spe] for spe in ancpopclade.get_leaf_labels()) - set(addStrtoTuple(ancpopname, ())))
	if verbose:
		print "ancpopname", ancpopname
		print "ancpopclade", repr(ancpopclade)
		print "lostpopnames", lostpopnames
	return (ancpopname, ancpopclade, lostpopnames)

def get_clade_phyloprofiles(dcolclades, spetree, **kw):
	ordspe = spetree.get_leaf_labels()
	dcolclaprof = {}
	for cla, leaflabs in dcolclades.iteritems():
		lspe = getSpeFromGenes(leaflabs, **kw)
		prof = [lspe.count(spe) for spe in ordspe]
		dcolclaprof[cla] = prof
	return dcolclaprof

def write_out_clade_phyloprofiles(dfamcolclaprof, nfout):
	fout = open(nfout, 'w')
	for fam in sorted(dfamcolclaprof.keys()):
		dcolclaprof = dfamcolclaprof[fam]
		for cla in sorted(dcolclaprof.keys()):
			fout.write('\t'.join([fam+'-'+cla]+[str(k) for k in dcolclaprof[cla]])+'\n')
	fout.close()

def getsavepoptree(nfpoptreeout="", poptree=None, spetree=None, lnamepops2collapse=None, speciestoprune=[], collapseAllPops=False, verbose=False):
	if not poptree: ptree = collapsePopulationInSpeciesTree(spetree, lnamepops2collapse, speciestoprune=speciestoprune, nested=True, collapseAllPops=collapseAllPops, verbose=verbose)
	else: ptree = poptree
	if nfpoptreeout:
		# save along the output gene tree sample file (species pruning being family-specific)
		cpoptree = copy.deepcopy(ptree)
		for node in cpoptree: # for compatibility with ALE
			if not node.is_leaf(): node.edit_label('')
		cpoptree.write_newick(nfpoptreeout)
	return ptree

def _allinterpopdist(tupi):
	i, lpopnames, poptree = tupi
	dnd = {}
	npop = len(lpopnames)
	popname1 = lpopnames[i]
	for j in range(i+1, npop):
		popname2 = lpopnames[j]
		nd = poptree[popname1].node_distance(poptree[popname2]) # min node dist is 2
		dnd[popname2] = nd
		if verbose: print popname1, popname2, nd
	return dnd

def allPWpopDist(poptree, lnamepops, nfmatpopdist, nbthreads=1, verbose=False):
	"""compute population all pairwise inter-node distances between populations on the (possibly collapsed) species tree"""
	if verbose: print 'allPWpopDist()'
	lpopnames = [n for n,p in lnamepops]
	npop = len(lpopnames)
	
	#~ def allinterpopdist(i):
		#~ dnd = {}
		#~ popname1 = lpopnames[i]
		#~ for j in range(i+1, npop):
			#~ popname2 = lpopnames[j]
			#~ nd = poptree[popname1].node_distance(poptree[popname2]) # min node dist is 2
			#~ dnd[popname2] = nd
			#~ if verbose: print popname1, popname2, nd
		#~ return dnd

	if nbthreads==1:
		ddnd = {}
		for i in range(npop):
			popname1 = lpopnames[i]
			ddnd[popname1] = _allinterpopdist((i, lpopnames, poptree))
	else:
		pool = mp.Pool(processes=nbthreads)
		res = pool.map(_allinterpopdist, ((i, lpopnames, poptree) for i in range(npop)))
		ddnd = {lpopnames[i]:dnd for i, dnd in enumerate(res)}
		
	with open(nfmatpopdist, 'w') as fmatpopdist:
		fmatpopdist.write('\t'.join(['']+lpopnames)+'\n')
		for popname1 in lpopnames:
			dnd = ddnd[popname1]
			fmatpopdist.write('\t'.join([popname1]+[str(dnd.get(popname2, ddnd[popname2].get(popname1, 0))) for popname2 in lpopnames])+'\n')
	return ddnd

def inferPopfromSpeTree(nfspetree, \
         populations=None, nfpopulationtree=None, nfdpopnodedist=None, \
         pop_stem_conds=None, within_pop_conds=None, \
         nbthreads=1, verbose=False, poptag=None):
	
	if verbose: print 'inferPopfromSpeTree()'
	# get file pre/sufix
	bnst, extst = nfspetree.rsplit('.', 1)
	## load species tree S
	spetree = tree2.read_check_newick(nfspetree, treeclass="AnnotatedNode")
	spetree.complete_internal_labels(order=0, ffel=True) # pre-order traversal labeling
	poptree = None
	## load pre-computed or compute populations in S
	if populations:
		if isinstance(populations, list):
			lnamepops = populations
		elif os.path.exists(populations):
			lnamepops = []
			if verbose: print 'load populations from \'%s\''% populations
			with open(populations, 'r') as fpop:
				for line in fpop:
					if line.startswith('#'): continue
					lsp = line.rstrip('\n').split('\t')
					lnamepops.append((lsp[0], tuple(lsp[1].split())))
		if nfpopulationtree:
			if verbose: print 'load population tree from \'%s\''% nfpopulationtree
			poptree = tree2.read_check_newick(nfpopulationtree, treeclass="AnnotatedNode")
			assert poptree.nb_leaves() == len(lnamepops)
	else:
		# define conditions to select species population within the species tree
		if pop_stem_conds: psc=pop_stem_conds
		else: psc = default_psc
		if within_pop_conds: wpc = within_pop_conds
		else: wpc = default_wpc
		# define species population that can be sued to summarize a collapsed clade
		lnamepops = speciesTreePopulations(spetree, pop_stem_conds=psc, within_pop_conds=wpc, poptag=poptag, testLeaves=True, nested=True, inclusive=False, taglen=3)
		# export in file ; name it after the species tree
		with open("%s_populations"%bnst, 'w') as foutspepop:
			foutspepop.write("# psc = %s\n"%repr(psc))
			foutspepop.write("# wpc = %s\n"%repr(wpc))
			for name, pop in lnamepops:
				foutspepop.write("%s\t%s\n"%(name, ' '.join(pop)))
	
	nfpoptreeout = "" if nfpopulationtree else "%s_collapsedPopulations.nwk"%(bnst)
	poptree = getsavepoptree(nfpoptreeout, poptree=poptree, spetree=spetree, lnamepops2collapse=lnamepops, collapseAllPops=True, verbose=verbose)
	poptree.prepare_fast_lookup()
	
	## annotate ancestor nodes of populations on S
	annotatePopulationInSpeciesTree(spetree, lnamepops)
	if nfpoptreeout:
		tnames, tpops = zip(*lnamepops)
		colour_tree_with_constrained_clades(spetree, tpops)
		spetree.write_nexus(nfpoptreeout.rsplit('.nwk', 1)[0]+'.nex', ignoreBS=True)
		tree2.dump_pickle(spetree, nfpoptreeout.rsplit('.nwk', 1)[0]+'.pickle')
		
	dspe2pop = getdspe2pop(lnamepops, spetree)
	
	## compute or load matrix of all pairwise population inter-node distances on the species tree
	if nfdpopnodedist:
		dpopnd = {}
		if verbose: print 'load pairwise population inter-node distances from \'%s\''% nfdpopnodedist
		with open(nfdpopnodedist, 'r') as fdpopnd:
			lpopnames = fdpopnd.readline().rstrip('\n').split('\t')[1:]
			for line in fdpopnd:
				lsp = line.rstrip('\n').split('\t')
				dpopnd[lsp[0]] = dict(zip(lpopnames, (int(d) for d in lsp[1:])))
		assert len(dpopnd) == len(lnamepops)
	else:
		dpopnd = allPWpopDist(poptree, lnamepops, "%s_interNodeDistPopulations"%(bnst), nbthreads=nbthreads, verbose=verbose)
	
	return (spetree, poptree, dspe2pop, lnamepops, dpopnd)
	
def mapPop2GeneTree(nfingt, dircons, dirout, method, spetree, poptree, dspe2pop, lnamepops, \
                    dpopnd={}, reuseOutput=0, combine_monophyletic=3, flatRCs=False, \
                    chain1ext='mb.run1.t', nbchains=2, aliext='nex', contreext='mb.con.tre', \
                    constrainttag='mbconstraints', collapsalntag='collapsed_alns', \
                    phyloproftag='phyloprofile', ccmaxlentag='max_subtree_lengths', \
					dirfullgt=None, dircolali=None, inalnext='nex', \
                    species_sep='_', charfilter=None, dontReplace=False, verbose=False, **kw):
	"""take as main input a name of Nexus tree chain file (from which is inferred 
	the list of all parallel chain files for the same family) and ouput a single 
	Newick tree chain file (ending '-Gtree.nwk') combining these parallel chains, 
	having replaced the collapsed clades with approapriate label or subtree.
	
	Depending on the method (see below), a tailored species tree will be generated
	in a file ending '-Stree.nwk', so that it covers the species set represented 
	in the output gene tree.
	
	Other files describing the data transformation are produced, including:
	- a phyloprofile (a presence/absence matrix of collapsed clade in species populations)
	- a table of correspondence (1-to-many mapping) between collapsed clades
	   and replacing leaf labels (file ending '-leaflabels_Spe2Pop.txt')
	- a table of what replacement approach led to the occurence of species in the ouput gene tree
	  (file ending '-replacedlabelsbyPop.txt')
	- a table of population occurence frequency in each colpased clades
	  (file ending '-PopFreqDistrib.txt')
	
	Using the input collapsed gene alignments found in the collpased alignment folder
	(in the collapsed gene tree info folder indicated by 'dircons', 
	the sub-folder with name indicated by 'collapsalntag'), this fonction
	will also generate alignments with sequences matching the replacement clades, 
	by duplicating the representative sequence so that the replacement clade 
	will be a flat rake in resulting gene tree.
		
	- if dontReplace is True, skip all the label/subtree replacement 
	  and just re-format the tree files from Nexus tree chain file(s)
	  to a single Newick tree chain file; only the 'Gtree.nwk' file is produced.
	
	- method can be 'collapseCCinG', 'replaceCCinGasinS', 'replaceCCinGasinS-collapsePOPinSnotinG'.
	
	With:
	 - P is the set of populations {Pi}, each being a set composed of 
	   member species leaves (pik) from the species tree S, mappping to 
	   a subtree si;
	 - C is the set of "bush" clades {Cj} in the gene tree G, defined by a
	   well-suported stem and no strong support for its topology within
	   (clades of very similar gene sequences with litlle or no available
	   phylogenetic information on their relationship);
	 - Pj = Pa(Cj) is the population identified as the ancestor of a "bush" clade Cj.
	 
	the proposed mehods are:
	
	- 'collapseCCinG': 
	will keep S as is, and replace all Cj in G by a single leaf labelled Pj
	# this option can be used to assign an ancestor species identity to gene leaf;
	# support for this is not yet implemented in reconciliation softwares like ALE or ecceTERA
		      
	- 'replaceCCinGasinS': 
	will keep S as is, and replace all Cj in G by a toppological equivalent 
	of the species subtree sj of the inferred ancestral population Pj, 
	and (if the original full gene tree is provided) scale the height of the 
	new gene subtree to that of the original Cj
		      
	- 'replaceCCinGasinS-collapsePOPinSnotinG': 
	same as above, but in addition collapse all Pi in S that are not 
	represented in the final form of G (once bush clades are replaced
	by single leaves or monophyletic species clades)
	
	- 'noreplace': do nothing; only to be used when no gene tree clade 
	has been collapsed in the first place.
	"""
	
	def replaceCCwithLabel(cla, speorpop, fam, dold2newname, tag='CC'):
		"""collapsed clade is replaced by a single leaf with a species or species population label
		
		new G leaf label has the following structure: POPNAME_FAM_XX-xxxx, 
		with XX=CC referring to a collapsed 'bush' clade; xxxx then contains the collapsed clade id e.g. 'clade1234'
		or XX=SG referring to a single gene sequence (located out of a collapsed 'bush' clade); xxxx then contains the original CDS code 'SPENAME_6789'.
		"""
		newname = "%s_%s_%s-%s"%(speorpop, fam, tag, cla)
		dold2newname[cla] = newname
		if verbose: print 'rename leaf:', newname
	
	def replaceCCwithSubtree(cla, ancpopnodelabels, fam, dold2newname, lostpops=[], speciestoprune=[], dccmaxlen={}, tag='RC', flatRCs=False, mode='tree2.Node'): #, onlyancs=False
		"""collapsed clade is replaced by a subtree of populations and/or single species copied from the population tree (collapsed species tree).
		
		new G subtree has leaf labels with the following structure: POPNAME_FAM_RC-xxxx, 
		POPNAME is the population or species label in the species tree; xxxx then contains the collapsed clade id e.g. 'clade1234'.
		"""
		if verbose: print "replaceCCwithSubtree(%s)"%repr(ancpopnodelabels)
		spesubtree = spetree[ancpopnodelabels] # ancpop can be a label (string) or a tuple of labels, in which case their MRCA will be fetched via Node.coalesce()
		sublnamepops = [(apnl, dnamepops.get(apnl, (apnl,))) for apnl in ancpopnodelabels]
		lostspecies = set(reduce(lambda x,y: x+y, (list(dnamepops[lp]) for lp in lostpops), []))
		spetoprune = list((set(speciestoprune) | lostspecies) & set(spesubtree.get_leaf_labels()))
		popsubtree = collapsePopulationInSpeciesTree(spesubtree, sublnamepops, speciestoprune=spetoprune, collapseAllPops=False, verbose=verbose)
		if verbose: print "replace with subtree of populations %s (%d leaves)"%(repr(ancpopnodelabels), popsubtree.nb_leaves())
		if flatRCs:
			# set the height of all branches in the new gene subtree to zero, i.e. a 'flat' replacement clade
			popsubtree *= 0
		elif dccmaxlen:
			# scale the height of the new gene subtree to that of the original
			maxlenspesub = popsubtree.max_leaf_distance()
			maxlenoricc = float(dccmaxlen[cla])
			if verbose: print 'maxlenspesub:', maxlenspesub, 'maxlenoricc', maxlenoricc
			if maxlenspesub>0:
				lenratio = maxlenoricc/maxlenspesub
				popsubtree *= lenratio
				if verbose: print 'multiplies branch length by factor %g'%lenratio
			elif maxlenoricc>0:
				# scale uniformly over the leaf branches
				if verbose: print 'set leaf branch to equal lengths %g'%maxlenoricc
				for popleaf in popsubtree:
					popleaf.set_lg(maxlenoricc)
		popsubtree.set_lg(0.0) # set its root length to 0
		for popnode in popsubtree: 
			if popnode.is_leaf():
				# new G leaf labels have the following structure: SPENAME_FAM_RC-clade1234, with RC referring to a replaced clade
				newspename = "%s_%s_%s-%s"%(popnode.label(), fam, tag, cla)
				popnode.edit_label(newspename)
			else:
				# naming internal nodes might lead to Gene being trees incompatible with ALE
				popnode.edit_label("")
		if verbose: print popsubtree.newick(ignoreBS=1)
		if mode=='tree2.Node':
			dold2newname[cla] = popsubtree
		elif mode=='Bio.Phylo':
			print "Warning: conversion of clade objects from tree2.Node to Bio.Phylo.TreeElement format is DEPRECATED in pantagruel/usingGeneRax branch"
			# already formats the tree in Bio.Phylo format to later integration to the parsed gene tree chain	
			dold2newname[cla] = tree2toBioPhylo(popsubtree)
		else:
			raise ValueError, "mode %s is not valid"%mode
	
	
	if verbose: print 'mapPop2GeneTree()'
	methods = ['noreplace', 'collapseCCinG', 'replaceCCinGasinS', 'replaceCCinGasinS-collapsePOPinSnotinG', \
	           'collapseCCinG-collapsematchingPOPinS', 'collapseALLinSandG', 'collapseCCinG-collapsePOPinSnotinGbutinCC']
	assert method in methods
	# get file pre/sufix
	dirgt = os.path.dirname(nfingt)
	bngt = os.path.basename(nfingt).split('.', 1)[0]
	outbn = bngt
	fam = bngt.rsplit('-', 1)[0]
	reptag = '' if dontReplace else '-replaced'
	
	nfoutcolGtree = "%s/%s%s-Gtree.nwk"%(dirout, outbn, reptag)
	if ((os.path.exists(nfoutcolGtree)) and (reuseOutput==2)):
		# assume work was done before, skip
		print "# (reused all) %s"%nfoutcolGtree
		return None
		
	if dontReplace:
		dold2newname = None
	else:
		dnamepops=dict(lnamepops)
		nfoutcolStree = "%s/%s-Stree.nwk"%(dirout, outbn)
		nfphyloprof = "%s_%s/%s.%s"%(dirout, phyloproftag, fam, phyloproftag)
		nfoutreflab = "%s/%s-leaflabels_Spe2Pop.txt"%(dirout, outbn)
		nfoutlabpoplab = "%s/%s-replacedlabelsbyPop.txt"%(dirout, outbn)
		nfoutfreqpopdistr = "%s/%s-PopFreqDistrib.txt"%(dirout, outbn)
		nfcolaln = "%s/%s/%s.%s"%(dircons, collapsalntag, outbn, inalnext)
		nfoutreplaln = "%s/%s-replaced.aln"%(dirout, outbn)
		nfccmaxlen = "%s/%s/%s.%s"%(dircons, ccmaxlentag, fam, ccmaxlentag)
		
		## define mapping from the leaf set
		ingt = tree2.read_newick(nfingt)
		gtleaflabels = ingt.get_leaf_labels()
		sspeG = set(getSpeFromGenes(gtleaflabels, species_sep=species_sep))
		genetreepopset = set([])
		for spe in sspeG:
			popname = dspe2pop.get(spe) # return None when original leaf label is not a canonical species, i.e. is a collapsed clade labelled 'cladeX'
			if popname:
				pop = dnamepops[popname]
				if len(pop)>1: genetreepopset.add(popname) # only consider non-trivial populations, i.e. grouping >= 2 species
		# record the set of populations absent in G;
		# will appear collapsed in S in methods '*-collapsePOPinSnotinG'
		notingttreepopset = (set(dnamepops.keys()) - genetreepopset)
		# load list of gene subtree constrained clades = clades collapsed in the full gene tree G
		consfilepat = '%s/%s/%s/%s*%s'%(dircons, constrainttag, fam, fam, constrainttag)
		if verbose: print "get collapsed clades from files matching pattern '%s'"%consfilepat
		lnfcons = glob.glob(consfilepat)
		if verbose: print repr_long_list(lnfcons)
		constraintclades = parseMrBayesConstraints(lnfcons)
		# write out the mapping of clades to species as phyloprofiles based on the species tree
		if ((not os.path.exists(nfphyloprof)) or (not reuseOutput)):
			famprof = get_clade_phyloprofiles(constraintclades, spetree)
			write_out_clade_phyloprofiles({fam:famprof}, nfphyloprof)
		# gather the max subtree lengths of the constrained clades
		dccmaxlen = {}
		if (constraintclades and os.path.exists(nfccmaxlen)) and (not flatRCs):
			if verbose: print "get max subtree lengths from file '%s'"%nfccmaxlen
			with open(nfccmaxlen, 'r') as fccmaxlen:
				for line in fccmaxlen:
					lsp = line.rstrip('\n').split('\t')
					dccmaxlen[lsp[0]] = float(lsp[1])
		elif dirfullgt:
			nffullgt = '%s/RAxML_rootedTree.%s.codes'%(dirfullgt, fam)
			fullgt = tree2.Node(file=nffullgt)
			# 'full' gene tree might actually not cover the whole set of leaves because RAxML was run on a non-redundant sequence alignment, i.e. with identical seq removed
			fullleafset = set(fullgt.get_leaf_labels()) 
			for cla, leaves in constraintclades.iteritems():
				nrleaves = fullleafset & set(leaves)
				if not nrleaves: raise ValueError, "no representative leaf of collapsed clade %s found in gene tree '%s'"%(cla, nffullgt)
				ccsullsubtree = fullgt.map_to_node(nrleaves)
				if (not flatRCs): dccmaxlen[cla] = ccsullsubtree.max_leaf_distance()
			if (not flatRCs):
				try:
					with open(nfccmaxlen, 'w') as fccmaxlen:
						fccmaxlen.write('\n'.join(['%s\t%f'%cml for cml in dccmaxlen.iteritems()])+'\n')
				except IOError, e:
					sys.stderr.write("failed attempt to write max subtree lengths of the constrained clades to file '%s'\n"%nfccmaxlen)

		# record collapsed bush clades' leaf labels to edit in G
		pop2collapseinS = set([])
		popnot2collapseinS = set([])
		dold2newname = {}
		speciestoprune = set([])
		foutfreqpopdistr = open(nfoutfreqpopdistr, 'w') if ((not os.path.exists(nfoutfreqpopdistr)) or (not reuseOutput)) else None
		for cla, leaves in constraintclades.iteritems():
			if verbose: print '##', cla
			lspeCC = getSpeFromGenes(leaves, species_sep=species_sep)
			dpopcount = _popCountDict(lspeCC, dspe2pop, dnamepops, asList=True)
			if foutfreqpopdistr: foutfreqpopdistr.write('\t'.join([cla]+['%s:%.3f'%(popname, mean(indcount)) for popname, indcount in dpopcount.iteritems()])+'\n')
			# will try and pool represented populations if they are monophyletic and have similar distributions of occurence
			if verbose: print "combineMonophyleticPopulations(%s)"%repr_long_list(dpopcount.keys())
			dpopcount = combineMonophyleticPopulations(dpopcount, poptree, dnamepops, dpopnd=dpopnd, maxpopnodedist=combine_monophyletic, **kw)
			# identify which species population is considered the ancestor of the collapsed clades
			if verbose: print "getCladeAncestor(%s)"%repr(lspeCC)
			ancpopname, ancpopclade, lostpopnames = getCladeAncestor(dpopcount, dnamepops, spetree, dspe2pop, verbose=verbose)
			# record the population(s) considered ancestral to the collasped clade
			# i.e. the one(s) used to replace the collasped clade in G (with a pop-labeled leaf or a S-tree-derived subtree) ;
			if type(ancpopname) is str:
				lancpopnames = [ancpopname]
			elif type(ancpopname) is tuple:
				lancpopnames = list(ancpopname)
			A = len(lancpopnames)
			pop2collapseinS |= set(lancpopnames) | set(lostpopnames)
			
			# collapsed clade is replaced by a S-tree-like subtree featuring the populations inferred as ancestrally present in the collapsed clade
			stleaves = set([])
			for t, pn in enumerate(lancpopnames+lostpopnames):
				if verbose: print ".%d %s"%(t, pn), ("(ancpop)" if t < A else "(lostpop)"),
				# evaluate presence of member species out of the clades to be colapsed
				ancpopmembers = set(dnamepops[pn])
				outCCancpopmembers = ancpopmembers & sspeG
				if verbose: print "; outCCancpopmembers:", outCCancpopmembers,
				if len(outCCancpopmembers)<=1:
					# if there are 0 or 1 members of the CC ancestral populations represented elsewhere in the Gene tree
					# better represent the CC with a single population label (CC tag)
					# if there is 1, outlier sequences of this member species will be equated with the population (SG tags, done below)
					if t < A:
						# this population is infered present in the ancestral clade, add to the subtree
						stleaves.add(pn)
				else:
					if t < A:
						# this population is infered present in the ancestral clade, add to the subtree
						# add a subtree covering the member species present out of the CC
						stleaves |= outCCancpopmembers
					# this population is not anymore to be collapsed in the Species tree
					popnot2collapseinS.add(pn)
					# record the remaining member species of the ancestral population for pruning in the species tree
					# i.e. those members NOT represented in the gene tree (even after substitution of CCs)
					speciestoprune = (speciestoprune | ancpopmembers) - sspeG
					if verbose: print "; speciestoprune:", speciestoprune,
				if verbose: print ""
			if len(stleaves)==1:
				replaceCCwithLabel(cla, stleaves.pop(), fam, dold2newname, tag='CC')
			elif len(stleaves)>1:
				replaceCCwithSubtree(cla, list(stleaves), fam, dold2newname, lostpops=lostpopnames, speciestoprune=list(speciestoprune), dccmaxlen=dccmaxlen, tag='RC')
			elif t < A:
				raise ValueError, "at least one species should be selected"
				
		if foutfreqpopdistr: foutfreqpopdistr.close()
		
		# define the set of populations of species absent in G, 
		# unless they are considered ancestral to any collapsed clade in G (including subsequent loss)
		# that will appear collapsed in S in methods '*-collapsePOPinSnotinG'
		pop2collapseinS = (pop2collapseinS | notingttreepopset) - popnot2collapseinS
		if verbose: print '#', fam, 'pop2collapseinS', pop2collapseinS
		
		if verbose: print '#', fam, 'dold2newname (CC)', dold2newname
		## prepare name reference dict for editing labels in gene tree sample
		dpop2SGnggen = {popname:0 for popname in dnamepops}
		if 'collapsePOPinS' in method:
			# collapse populations in S
			lnamepops2collapse = [(name, pop) for name, pop in lnamepops if (name in pop2collapseinS)]
		# S tree with some populations collapsed as ancestors and some not, used later for reconciliation inference on this particular family
		if verbose: print '#', fam, 'lnamepops2collapse:\n %s'%repr_long_list(lnamepops2collapse)
		# add species from populations not to collapse to the list of populations as standalone species
		dnamepops2collapse = dict(lnamepops2collapse)
		for name, pop in lnamepops:
			if name not in dnamepops2collapse:
				for spe in pop:
					if spe not in speciestoprune:
						lnamepops2collapse.append((spe, (spe,)))
							
		getsavepoptree(nfoutcolStree, spetree=spetree, lnamepops2collapse=lnamepops2collapse, speciestoprune=list(speciestoprune), verbose=max(0, verbose-1))
		# then record other single leaf labels belonging to collapsed pops in S to edit in G
		# (species that are the single representatives of collapsed populations in G)
		dspe2pop2collapse = getdspe2pop(lnamepops2collapse)
		for leaflab in gtleaflabels:
			spelab = leaflab.split(species_sep, 1)[0]
			poplab = dspe2pop2collapse.get(spelab)
			if poplab and poplab!=spelab:
				replaceCCwithLabel(leaflab, poplab, fam, dold2newname, tag='SG')
				dpop2SGnggen[poplab] += 1
					
		reptypes = ['CC', 'RC', 'SG']
		if verbose: print '#', fam, 'dold2newname (%s)'%('+'.join(reptypes)), dold2newname
		## output label translation reference
		dpoptyperepcount = {}
		if ((not os.path.exists(nfoutreflab)) or (not reuseOutput)):
			with open(nfoutreflab, 'w') as foutreflab:
				for oldlab, newlaborst in dold2newname.iteritems():
					newlabs = labsFromReplacementLabOrSubtree(newlaborst)
					for newlab in newlabs:
						newlabpop, fam, ngtag = newlab.rsplit('-'+oldlab, 1)[0].rsplit('_',2)
						dpoptyperepcount.setdefault(newlabpop, {rt:0 for rt in reptypes})[ngtag] += 1
						foutreflab.write('\t'.join((oldlab, newlab))+'\n')
		if dpoptyperepcount or ((not os.path.exists(nfoutlabpoplab)) or (not reuseOutput)):
			with open(nfoutlabpoplab, 'w') as foutlabpoplab:
				foutlabpoplab.write('\t'.join(['Pop_or_Spe_id']+reptypes)+'\n')
				for pop in dpoptyperepcount:
					foutlabpoplab.write('\t'.join([pop]+[str(dpoptyperepcount[pop][reptype]) for reptype in reptypes])+'\n')
			
	## apply to gene tree sample
	# load gene tree sample from Newick tree file, replace labels/clades in it, and write it out
	replaceInSingleTree(nfingt, dold2newname, nfgtout=nfoutcolGtree, verbose=verbose, maskchars=charfilter)
	# generate alignment with sequence labels matching the ones in the gene tree produced above;
	# for that, the sequece representative of the CC in the input collapsed alignment is duplicated 
	# into identical sequences with different labels as in the replacement clade 
	duplicateSeqsInAln(nfcolaln, dold2newname, nfoutreplaln, verbose=verbose)
		
def usage():
	s =  'For population assignation:\n'
	s += ' python2.7 replace_species_by_pop_in_gene_trees.py -S /path/to/species_tree\n'
	s += 'For replacement of species labels with population labels in gene tree samples:\n'
	s += ' python2.7 replace_species_by_pop_in_gene_trees.py -S /path/to/[TIMED.]species.tree -G /path/to/list.of.gene.tree.files -c /path/to/folder.of.collapsed.clade.info.files -o /path/to/output.folder \\\n'
	s += ' [--populations=/path/to/species2population_map --population_tree=/path/to/population_annotated_species_tree] [OTHER OPTIONS]\n'
	s += 'For NO replacement of species labels in gene tree samples (simply reformating):\n'
	s += ' python replace_species_by_pop_in_gene_trees.py -G /path/to/list.of.gene.tree.files --no_replace -o /path/to/output.folder\n'
	s += ''
	s += ' Options:\n'
	s += '  Input file options:\n'
	s += '  --populations\t\t\tpath to file with already computed mapping from species populations labels to species leaf labels.\n'
	s += '  --population_tree\t\t\tpath to Newick-formatted tree file of an already computed mapping of populations ancestors onto the species tree.\n'
	s += '  --population_node_distance\t\tpath to pre-computed matrix of inter-node distances between populations in the species tree.\n'
	s += '  --dir_full_gene_trees\t\tpath to folder of original (before clade collapsing) gene trees (can be used to scale subtree grafts)\n\t\t\t\t(overridden by the presence of a adequate file in the collapsed.clade.info.files subfolder \'max_subtree_lengths\').\n'
	s += '  Parameters:\n'
	s += '  --pop_stem_conds\t\tconditions on clade stem to select populations in the species tree (see mark_unresolved_clades.py)\n'
	s += '  --within_pop_conds\t\tconditions on clade subtree to select populations in the species tree (see mark_unresolved_clades.py)\n'
	s += '  --pop_lab_prefix\t\t(preferably short) string that will be systematically used as prefix for population names\n'
	s += '  --'+mapPop2GeneTree.__doc__+'\n'
	s += '  --threads\t\t\tnumber of parralel processes to run.\n'
	s += '  --reuse\t\t[0,1,2]\tif 0: normal behaviour;\n'
	s += '	       \t\t\t\tif 1: skip computing and writing output files which name matches already existing file;\n'
	s += '	       \t\t\t\t      BEWARE! This can save a lot of time but may produce inconsistent output file sets!\n'
	s += '	       \t\t\t\tif 2: skip processing the entire gene tree chain if the expected output gene tree chain file (with labels replaced) exists;\n'
	s += '  --max.recursion.limit\t(int)Set maximum recursion limit of Python, may be required for processing big trees  (default 4000).\n'
	s += '  --no_replace won\'t touch the tree structure, just reformat them from several Nexus tree chain files into a single Newick chain tree file.\n'
	s += '  --logfile\t\t\tpath to file where combined stdout and stderr are redirected (if multiple threads, several files with this basename and a pid suffix will be created in addition to the master process\' one).\n'
	s += '  --filter_dashes\t\t\tenables filtering of dashes PRECEDED BY A CAPITAL LETTER OR DIGIT when parsing the Nexus trees; dash are restored when written to Newick.\n'
	s += '  --verbose {0,1,2}\tverbose mode, from none to plenty.\n'
	s += '  -v\t\t\tequivalent to --verbose=1.\n'
	return s

	
if __name__=='__main__':	
	opts, args = getopt.getopt(sys.argv[1:], 'G:c:S:o:T:hv', ['dir_out=', 'method=', 'no_replace', 'flat_RCs', \
															'populations=', 'population_tree=', 'population_node_distance=', 'dir_full_gene_trees=', \
															'pop_stem_conds=', 'within_pop_conds=', 'pop_lab_prefix=', 'filter_dashes', \
															'threads=', 'reuse=', 'verbose=', 'max.recursion.limit=', 'logfile=', 'help'])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	nfspetree = dopt.get('-S')
	nflnfingt = dopt.get('-G')
	dircons = dopt.get('-c')
	dirout = dopt.get('-o')
	dontReplace = ('--no_replace' in dopt)
	defmethods = 'replaceCCinGasinS-collapsePOPinSnotinG' if (not dontReplace) else 'noreplace'
	methods = dopt.get('--method', defmethods).split(',')
	nfpops = dopt.get('--populations')
	nfpoptree = dopt.get('--population_tree')
	nfdpopnd = dopt.get('--population_node_distance')
	dirfullgt = dopt.get('--dir_full_gene_trees')
	psc = eval(dopt.get('--pop_stem_conds', 'None'))
	wpc = eval(dopt.get('--within_pop_conds', 'None'))
	reuseOutput = int(dopt.get('--reuse', 0))
	verbose = int(dopt.get('--verbose', ('-v' in dopt)))
	maxreclim = int(dopt.get('--max.recursion.limit', 4000))
	nflog = dopt.get('--logfile')
	nbthreads = int(dopt.get('--threads', dopt.get('-T', -1)))
	poptag = dopt.get('--pop_lab_prefix')
	flatRCs = ('--flat_RCs' in dopt)
	if ('--filter_dashes' in dopt):
		charfilter = [('\([A-Z0-9]\)-', '\\1@'), ('\([A-Z0-9]\)@', '\\1-')]
	else:
		charfilter = None
	if nbthreads < 1: nbthreads = mp.cpu_count()
	
	sys.setrecursionlimit(maxreclim)
	
	if nflog:
		flog = open(nflog, 'w')
		sys.stdout = sys.stderr = flog
	
	if not dontReplace:
		if not nfspetree:
			raise ValueError, "Providing a species tree file is mandatory when replacing clades in gene tree chains"
		### define species populations
		spetree, poptree, dspe2pop, lnamepops, dpopnd = inferPopfromSpeTree(nfspetree, populations=nfpops, poptag=poptag, \
                                                               nfpopulationtree=nfpoptree, nfdpopnodedist=nfdpopnd, \
                                                               pop_stem_conds=psc, within_pop_conds=wpc, \
                                                               nbthreads=nbthreads, verbose=verbose)
	else:
		spetree = poptree = dspe2pop = lnamepops = dpopnd = None
	
	if nflnfingt:
		ldircreate = []
		if not dirout: raise ValueError, "must provide value for options -o when passing gene tree samples to -G"
		ldircreate += [os.path.join(dirout, meth) for meth in methods]
		if not dontReplace:
			if not dircons: raise ValueError, "must provide value for options -c when passing gene tree samples to -G for label replacement (i.e. without option --no_replace)"
			ldircreate += ['_'.join([os.path.join(dirout, meth), phyloproftag]) for meth in methods]
		else:
			dircons = None
		for d in ldircreate:
			if not os.path.isdir(d):
				os.mkdir(d)
		
		### map populations on gene trees	
		with open(nflnfingt, 'r') as flnfingt:
			lnfingt = [line.strip('\n') for line in flnfingt]
		
		# wrapper function for a single-family run
		def mapPop2GeneTreeSetArgs(nfingt):
			if nflog and nbthreads>1:
				pflog = open(nflog+".%d"%(os.getpid()), 'a')
				sys.stdout = sys.stderr = pflog
				
			print '#', nfingt
			for method in methods:
				print 'method:', method
				try:
					mapPop2GeneTree(nfingt, dircons, os.path.join(dirout, method), method, spetree, poptree, dspe2pop, lnamepops, dpopnd, \
									dontReplace=dontReplace, reuseOutput=reuseOutput, flatRCs=flatRCs, phyloproftag=phyloproftag, dirfullgt=dirfullgt, charfilter=charfilter, verbose=verbose)
				except Exception, e:
					if (type(e) is RuntimeError) and str(e).startswith('maximum recursion depth exceeded'):
						print "met RuntimeError on file '%s':"%nfingt, e
						if sys.getrecursionlimit() < maxreclim*maxreclimupfactor:
							# allow for increase of recursion limit
							sys.setrecursionlimit(sys.getrecursionlimit()+2000)
							print "re-run with increased recursion depth +2000"
							mapPop2GeneTreeSetArgs(nfingt)
						else:
							print "Error: failed to run on file '%s':"%nfingt
							return 1
					else:	
						print 'Caught exception:'
						traceback.print_exc()
						sys.stdout.flush()
						raise e
			
			print '# # # #\n'
			if nflog and nbthreads>1:
				if not pflog.closed:
					pflog.close()
			return 0
		
		# run in sequential
		if nbthreads==1:
			for nfingt in lnfingt:
				if verbose: print nfingt
				mapPop2GeneTreeSetArgs(nfingt)
		# or run in parallel
		else:
			pool = mp.Pool(processes=nbthreads)
			res = pool.map(mapPop2GeneTreeSetArgs, lnfingt)
			#~ res = [pool.apply_async(mapPop2GeneTreeSetArgs, nfingt) for nfingt in lnfingt] #, chunksize=5

	if nflog:
		flog.close()
