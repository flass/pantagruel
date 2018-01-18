#!/usr/bin/python

import tree2
import sys, os, getopt, glob
from mark_unresolved_clades import select_clades_on_conditions, mean, median, var
from scipy import stats
import copy

import multiprocessing as mp
from Bio import File, AlignIO, Align
from Bio.Phylo import BaseTree, NewickIO, NexusIO
from Bio.Phylo import _io as PhyloIO
from StringIO import StringIO
import traceback

sys.setrecursionlimit(4000)

supported_formats = {'newick': NewickIO, 'nexus': NexusIO}
phyloproftag='phyloprofiles'

def repr_long_list(l, lenbeg=5, lenend=5):
	assert type(l) is list
	if len(l) <= lenbeg+lenend:
		return repr(l)
	else:
		if lenbeg>0: a = '%s, '%( repr(l[:lenbeg]).strip('[]') )
		else: a = ''
		if lenend>0: b =', %s'%( repr(l[(-1*lenend):]).strip('[]') )
		else: b = ''
		return '[%s...%s] (%d items)'%(a, b, len(l))

def parseMrBayesConstraints(lnfcons):
	constraintclades = {}
	for nfcons in lnfcons:
		with open(nfcons, 'r') as constraintsfile:
			for line in constraintsfile:
				if line.startswith('constraint'):
					scla, sep, sleaves = line.rstrip(';\n').partition(' = ')
					cla = scla.split()[1]
					leaves = sleaves.split()
					constraintclades[cla] = leaves
	return constraintclades

def treeappend(trees, file, format, **kwargs):
    """appending version of Phylo._io.write()"""
    if isinstance(trees, (BaseTree.Tree, BaseTree.Clade)):
        # Passed a single tree instead of an iterable -- that's OK
        trees = [trees]
    with File.as_handle(file, 'a+') as fp:
        n = getattr(supported_formats[format], 'write')(trees, fp, **kwargs)
    return n

def in2outChainFileName(nfchain, dirout, inchainfmt='nexus', outchainfmt='newick', **kw):
	dfmtext = {'nexus':'nex', 'newick':'nwk'}
	if os.path.isdir(dirout):
		nfout = "%s/%s"%(dirout, dfmtext[outchainfmt].join(os.path.basename(nfchain).split(dfmtext[inchainfmt])))
	else:
		nfout = "%s.%s"%(dirout, dfmtext[outchainfmt])
	return nfout

def tree2toBioPhylo(node):
	"""perform tree2.Node -> Newick -> Bio.Phylo.Newick.Tree conversion"""
	return PhyloIO.read(StringIO(node.newick()), 'newick', rooted=True)

def parseChain(lnfchains, dold2newname, nfchainout=None, inchainfmt='nexus', outchainfmt='newick', verbose=False):
	"""parse a Nexus-format tree chain and re-write it in a Newick format; edit the tree on the fly, substituting tip labels or grafting trees on tips
	
	
	"""
	if verbose: print "parseChain('%s')"%repr(lnfchains)
	ntree = 0
	buffsize = 50
	treebuffer = []
	if nfchainout: nfout = nfchainout
	else: nfout = n2outChainFileName(lnfchains[0], dirout, inchainfmt, outchainfmt)
	dhandles = {}
	for k, nfchain in enumerate(lnfchains):
		dhandles[k] = PhyloIO.parse(nfchain, inchainfmt)
	# intertwines the trees from separate chains, respecting their sequence
	# so that the burn-in can still be defined as the first portion of trees in the sample.
	sometreeleft = True
	ichains = range(len(lnfchains))
	while sometreeleft:
		for k in ichains:
			#~ sys.stdout.write('\rk'+str(k)) ; sys.stdout.flush()
			try:
				tree = dhandles[k].next()
				#~ print tree
				#~ sys.stdout.flush()
			except StopIteration:
				sometreeleft = False
				# stops for all chains here, even if they had different length (which shoulld not be)
				break # the for k loop
			if dold2newname:
				for tip in tree.get_terminals():
					nutipname = dold2newname.get(tip.name)
					if nutipname:
						if isinstance(nutipname, str):
							tip.name = nutipname
						else:
							if isinstance(nutipname, tree2.Node):
								nusubtree = tree2toBioPhylo(nutipname)
							elif isinstance(nutipname, BaseTree.TreeElement):
								nusubtree = nutipname
							else:
								raise ValueError, "replacement value for a leaf must be either a string (to edit leaf label) of a tree object instance of classes tree2.Node or Bio.Phylo.TreeElement (or derivates)"
							subtreelen = nusubtree.total_branch_length()/nusubtree.count_terminals()
							if subtreelen:
								# substract the subtree length to its branch length
								tip.branch_length = max(0.0, tip.branch_length - subtreelen)
							# attach subtree to tree
							tip.clades = nusubtree.root.clades
							tip.name = ''
			treebuffer.append(tree)
			ntree += 1
			#~ sys.stdout.write('\rntree'+str(ntree)) ; sys.stdout.flush()
			if len(treebuffer) >= buffsize:
				if verbose:
					sys.stdout.write('\r'+str(ntree)) ; sys.stdout.flush()
				if ntree <= buffsize:
					PhyloIO.write(treebuffer, nfout, outchainfmt)
				else:
					treeappend(treebuffer, nfout, outchainfmt)
				treebuffer = []
	if treebuffer: treeappend(treebuffer, nfout, outchainfmt)
	print '%s%s ...done (%d trees)'%(('\n' if verbose else ''), nfout, ntree)
	return None

def getdspe2pop(lnamepops, spetree=None):
	dspe2pop = {}
	for popname, pop in lnamepops:
		for spe in pop:
			dspe2pop[spe] = popname
	if spetree:
		# completes with species on their own in species tree
		for spe in spetree.get_leaf_labels():
			if spe not in dspe2pop:
				dspe2pop[spe] = spe
				lnamepops.append((spe, (spe,)))
	return dspe2pop

def annotatePopulationInSpeciesTree(spetree, lnamepops, returnCopy=False, returnAncNodes=False):
	"""use when species populations and their single species members can coexist in the species/gene trees.
	
	i.e. if reconciliation can handle a gene tree with both extant  and ancestral species at its leaves.
	"""
	if returnCopy: poptree = copy.deepcopy(spetree)
	else: poptree = spetree
	lanc = []
	for popname, pop in lnamepops:
		popanc = poptree.mrca(pop)
		popanc.edit_label(popname)
		lanc.append(popanc)
	if returnCopy: return poptree
	elif returnAncNodes: return lanc
	else: return None

def speciesTreePopulations(spetree, pop_stem_conds, within_pop_conds, nested=False, inclusive=False, taglen=3, **kw):
	"""find supported clades in the species tree
	
	when inclusive is False, populations can be defined so as to emerges from another, with the more ancestral not be inclusive on of another, i.e. a population
	"""
	#~ lpops = mark_unresolved_clades(spetree , clade_stem_conds=pop_stem_conds, within_clade_conds=within_pop_conds, threshold_method='fixed', nested=nested, inclusive=inclusive, **kw)
	lpops = select_clades_on_conditions(spetree , clade_stem_conds=pop_stem_conds, within_clade_conds=within_pop_conds, testRoot=True, nested=nested, inclusive=inclusive, **kw)
	# generate pop names based on the major 3 first letter of species in it
	lnamepops = []
	dnamecount = {}
	for pop in lpops:
		if len(pop)>1:
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

def collapsePopulationInSpeciesTree(spetree, lnamepops, fast=True, speciestoprune=[], nested=True, keepUltrametric=True, collapseAllPops=False, verbose=0):
	"""if population are nested, the order of the input population (name, (members,)) list matters!
	verbose=0: silent, verbose=1: general progression, verbose=2: detail of tree edition.
	"""
	if verbose: print 'collapsePopulationInSpeciesTree()'
	# first mark the population ancestor node
	spetree.complete_node_ids(order=2) # post-order traversal numeration
	poptree = copy.deepcopy(spetree)
	dnamepop = dict(lnamepops)
	# list of nodes independent of tree structure (i.e. pruning or unlinking nodes in the tree won't affect list sequence]
	lnamepopancs = [(popname, poptree.mrca(popspe)) for popname, popspe in lnamepops]
	if verbose: print 'lnamepopancs:\n%s'%( repr_long_list(lnamepopancs) )
	if fast:
		# ensure sorting pop. ancestor nodes following post-order traversal numeration in the original tree
		lnamepopancs.sort(key=lambda x: x[1].nodeid())
		if verbose: print 'lnamepopancs (sorted according to popancs nodeids):\n%s'%( repr_long_list(lnamepopancs) )
	lpopname, lpopancs = zip(*lnamepopancs)
	lpopancids = [popanc.nodeid() for popanc in lpopancs] # increasing sequence of node ids
	for a, (popname, popanc) in enumerate(lnamepopancs):
		lpopspe = dnamepop[popname]
		if verbose: print '#%d. %s:'%(a, popname), lpopspe
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
				# <=> [k:a] with a = lpopancids.index(p) and k = min([t for t, j in enumerate(lpopancids) if j>=(a-n+1)])
				n = spetree.idgetnode(popanc.nodeid()).nb_all_children() # rely on parent-child structure of the original, uncut tree
				k = 0
				for k, j in enumerate(lpopancids):
					if j>=(a-n+1): break
				if verbose: print 'check subancs within lpopancs[%d:%d]'%(k,a)
				subancs = popanc.is_parent_of_any(lpopancs[k:a], returnList=True)
			else:
				if verbose: print 'check all subancs'
				subancs = popanc.is_parent_of_any(lpopancs, returnList=True)
			if verbose: print 'subancs (labels):', [su.label() for su in subancs]
			if subancs:
				# remove population clade (or paraphyletic group) leaf by leaf
				# until only other population ancestor and one arbitrary leaf are left
				nkept = 0
				for spe in lpopspe:
					if spe in popanc.children_labels():
						# leaf is just below the population ancestor 
						# must not pop() that leaf to avoid destroying the population ancestor node
						if nkept==0:
							# meeting the population ancestor must happen at least once, for which the leaf is kept
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
						# prune the leaf (using tree2.Node.pop()):
						# lead to deconnecting its father node from the rest of tree 
						# and connecting the brother node to the grand-father node
						if verbose: print 'pop', spe #,
						pspe = popanc.pop(spe)
						#~ print pspe
	
			else:
				# no other population ancestor below
				#~ assert set(popanc.get_leaf_labels())==set(lpopspe)
				if set(popanc.get_leaf_labels())!=set(lpopspe):
					print 'set(popanc.get_leaf_labels()):'
					print set(popanc.get_leaf_labels())
					print 'set(lpopspe):'
					print set(lpopspe)
					print 'set(popanc.get_leaf_labels()) - set(lpopspe)'
					print set(popanc.get_leaf_labels()) - set(lpopspe)
					raise ValueError, "there should be no other population ancestor below"
				popanc.as_leaf(newlabel=popname, silent=1-(verbose-1))
				if keepUltrametric and sublen!='NA': popanc += sublen	# reports distances to longest leaf lineage below the ancestor so that an ultrametric tree keeps being ultrametric
				if verbose: print 'as_leaf', popname, '+=', sublen
						
		else:
			# assume populations are monophyletic and the clade can be removed at once
			assert set(popanc.get_leaf_labels())==set(lpopspe)
			popanc.as_leaf(newlabel=popname, silent=1-(verbose-1))

	if collapseAllPops and (not set(lpopname) == set(poptree.get_leaf_labels())):
		print "set(lpopname) - set(poptree.get_leaf_labels()):\n", set(lpopname) - set(poptree.get_leaf_labels())
		print "set(poptree.get_leaf_labels()) - set(lpopname):\n", set(poptree.get_leaf_labels()) - set(lpopname)
		raise IndexError, "population list and population tree are discordant"
		
	if speciestoprune:
		for spetoprune in speciestoprune:
			s2p = poptree[spetoprune]
			if s2p: poptree.pop(s2p)
			
	return poptree

#~ nfspetree = '/home/flassall/PanteroDB/core_genome/raxml_tree/RAxML_bipartitions.concat_core_cds_464entero.rooted.codes.nwk'
#~ spetree = tree2.AnnotatedNode(file=nfspetree)
#~ spetree.complete_internal_labels(order=0, ffel=True) # pre-order traversal labeling
#~ lnamepops = speciesTreePopulations(spetree, pop_stem_conds=[('lg', '>=', 0.0002), ('bs', '>=', 80)], within_pop_conds=[('max', 'lg', '<=', 0.0002, -1)], testLeaves=True, nested=True, inclusive=False, taglen=3)
#~ dspe2pop = getdspe2pop(lnamepops, spetree)
#~ poptree = collapsePopulationInSpeciesTree(spetree, lnamepops, nested=True)
#~ poptree.figtree()
#~ annotatePopulationInSpeciesTree(spetree, lnamepops)

def _popCountDict(lspe, dspe2pop, dnamepops, asList=True):
	dspecount = {spe:lspe.count(spe) for spe in set(lspe)}
	dpopcount = {}
	if asList:
		# will return occurence profile (count per strain)
		popnames = set([dspe2pop[spe] for spe in dspecount])
		for popname in popnames:
			dpopcount[popname] = [dspecount.get(spe, 0) for spe in dnamepops[popname]]
	else:
		# will return total occurence count
		for spe in dspecount:
			popname = dspe2pop[spe]
			dpopcount.setdefault(popname, 0) + dspecount[spe]
	return dpopcount

def getSpeFromGenes(leaflabs, species_sep='_', **kw):
	return [leaflab.split(species_sep)[0] for leaflab in leaflabs]

def getCladeAncestor(lspe, spetree, dspe2pop, dnamepops, method='max_meanfreqpop', pvaltresh=0.05, dpopnd={}, **kw):
	verbose = kw.get('verbose', False)
	## determine who was there first
	crit, estim = method.split('_')
	## based on a census of represented species
	dpopcount = _popCountDict(lspe, dspe2pop, dnamepops, asList=True)
	if estim=='totalcountpop': foo = sum
	elif estim=='meanfreqpop': foo = mean
	elif estim=='medianfreqpop': foo = median
	else: raise ValueError, "chose wrong method for ranking puative ancestor populations"
	allpopdistr = [(popname, float(foo(popcount)), popcount) for popname, popcount in dpopcount.iteritems()]
	if len(allpopdistr) > 1:
		## first rank populations according to estimator of targeted quantity and given optimal criterion
		if crit=='max':
			a = -1
		elif crit=='min':
			a = 1
		else:
			raise ValueError
		# always from best values (index 0) to worst (index -1), breaking ties with sample size (more = better)
		allpopdistr.sort(key=lambda x: (a*x[1], -1*len(x[2]))) 
		# get min/max population
		popname0, popestim0, indcounts0 = allpopdistr[0]
		# explore populations ranked second, etc. to test difference of distributions
		for k in range(1, len(allpopdistr)):
			popnamek, popestimk, indcountsk = allpopdistr[k]
			if verbose: print indcounts0, 'vs', indcountsk
			tt = stats.ttest_ind(indcounts0, indcountsk, equal_var=False)
			if verbose: print "difference of mean frequencies in %s (%f, n=%d) and %s (%f, n=%d): %s"%(popname0, popestim0, len(indcounts0), popnamek, popestimk, len(indcountsk), str(tt))
			if tt[1] < pvaltresh:
				# statistically different 
				break
		# all the first k populations have statistically equivalent frequencies
		## will chose the most central population in the species tree
		eqpopdistr = allpopdistr[:k]
		if len(eqpopdistr) > 1:
			if not dpopnd:
				dnd = {}
				# compute all pairwise node distance (on the species tree) between the populations
				for i, epd1 in enumerate(eqpopdistr):
					p1 = epd1[0]
					for j, epd2 in enumerate(eqpopdistr[i+1:]):
						p2 = epd2[0]
						nd = spetree[p1].node_distance(spetree[p2]) # min node dist is 2
						dnd.setdefault(p1, []).append(nd) 
						dnd.setdefault(p2, []).append(nd) 
				rankedpopdist = sorted([(ipop, sum(dnd[epd[0]])) for ipop, epd in enumerate(eqpopdistr)], key=lambda x: x[1])
			else:
				rankedpopdist = sorted([(ipop, sum([dpopnd[epd1[0]][epd2[0]] for epd2 in eqpopdistr])) for ipop, epd1 in enumerate(eqpopdistr)], key=lambda x: x[1])
			# find the population the least distant from the others
			
			ancpopname, ancpopestim, indcounts = eqpopdistr[rankedpopdist[0][0]]
			if verbose: print rankedpopdist, "most central population in tree is:", ancpopname
			
		else:
			ancpopname, ancpopestim, indcounts = eqpopdistr[0]
			
	else:
		ancpopname, ancpopestim, indcounts = allpopdistr[0]
	if estim.endswith('freqpop'):
		# keep track of roughly how many gene copies were there in the clade
		# to avoid biasing to much the reconciliation if there was 
		# a conserved ancestral duplication above the clade's node
		ancnbcopies = ancpopestim
	else:
		ancnbcopies = 1
	return (ancpopname, ancnbcopies, allpopdistr)

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

def getsavepoptree(nfpoptreeout, poptree=None, spetree=None, lnamepops2collapse=None, speciestoprune=[], collapseAllPops=False, verbose=False):
	if not poptree: poptree = collapsePopulationInSpeciesTree(spetree, lnamepops2collapse, speciestoprune=speciestoprune, nested=True, collapseAllPops=collapseAllPops, verbose=verbose)
	# save along the output gene tree sample file (species pruning being family-specific)
	cpoptree = copy.deepcopy(poptree)
	for node in cpoptree: # for compatibility with ALE
		if not node.is_leaf(): node.edit_label('')
	cpoptree.write_newick(nfpoptreeout)
	return poptree

def allPWpopDist(poptree, lnamepops, nfmatpopdist, nbthreads=1, verbose=False):
	"""compute population all pairwise inter-node distances between populations on the (possibly collapsed) species tree"""
	if verbose: print 'allPWpopDist()'
	lpopnames = [n for n,p in lnamepops]
	npop = len(lpopnames)
	
	def allinterpopdist(i):
		dnd = {}
		popname1 = lpopnames[i]
		for j in range(i+1, npop):
			popname2 = lpopnames[j]
			nd = poptree[popname1].node_distance(poptree[popname2]) # min node dist is 2
			dnd[popname2] = nd
			if verbose: print popname1, popname2, nd
		return dnd

	if nbthreads==1:
		ddnd = {}
		for i in range(npop):
			popname1 = lpopnames[i]
			ddnd[popname1] = allinterpopdist(i)
	else:
		pool = mp.Pool(processes=nbthreads)
		res = pool.map(allinterpopdist, range(npop))
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
         nbthreads=1, verbose=False):
	
	if verbose: print 'inferPopfromSpeTree()'
	# get file pre/sufix
	bnst, extst = nfspetree.rsplit('.', 1)
	## load species tree S
	spetree = tree2.AnnotatedNode(file=nfspetree)
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
			poptree = tree2.AnnotatedNode(file=nfpopulationtree)
			assert poptree.nb_leaves() == len(lnamepops)
	else:
		# define conditions to select species population within the species tree
		if pop_stem_conds: psc=pop_stem_conds
		else: psc = [('lg', '>=', 0.0002), ('bs', '>=', 80)]
		if within_pop_conds: wpc = within_pop_conds
		else: wpc = [('max', 'lg', '<=', 0.0002, -1)]
		# define species population that can be sued to summarize a collapsed clade
		lnamepops = speciesTreePopulations(spetree, pop_stem_conds=psc, within_pop_conds=wpc, testLeaves=True, nested=True, inclusive=False, taglen=3)
		# export in file ; name it after the species tree
		with open("%s_populations"%bnst, 'w') as foutspepop:
			foutspepop.write("# psc = %s\n"%repr(psc))
			foutspepop.write("# wpc = %s\n"%repr(wpc))
			for name, pop in lnamepops:
				foutspepop.write("%s\t%s\n"%(name, ' '.join(pop)))
				
	poptree = getsavepoptree("%s_collapsedPopulations.nwk"%(bnst), poptree=poptree, spetree=spetree, lnamepops2collapse=lnamepops, collapseAllPops=True, verbose=verbose)
	
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
	
	## annotate ancestor nodes of populations on S
	annotatePopulationInSpeciesTree(spetree, lnamepops)
	dspe2pop = getdspe2pop(lnamepops, spetree)
	
	return (spetree, dspe2pop, lnamepops, dpopnd)
	
def mapPop2GeneTree(nfingtchain1, dircons, dirout, method, spetree, dspe2pop, lnamepops, dpopnd={}, reuseOutput=0, \
                    chain1ext='mb.run1.t', nbchains=2, aliext='nex', contreext='mb.con.tre', \
                    constrainttag='mbconstraints', collapsalntag='collapsed_alns', phyloproftag='phyloprofile', ccmaxlentag='max_subtree_lengths', dirfullgt=None, \
                    species_sep='_', verbose=False, **kw):
	"""method can be 'collapseCCinG', 'replaceCCinGasinS', 'replaceCCinGasinS-collapsePOPinSnotinG, 'collapseinSandG', 'collapseALLinSandG' or 'collapseALLinSbutinGbutinCC.
	
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
	
	#                !!! the following methods are DEPRECATED !!!
	# they lead to change of G-to-S mapping cardinality, 
	# i.e. some leaves representing different species assigned to a same population 
	# but paraphyletic in the gene tree, or even monophyletic but not grouped under a same bush calde
	# will artefactually be seen as multiple gene copies for one single species population
	     
	- 'collapseCCinG-collapsematchingPOPinS': 
	for all Cj, collapse both Cj in G and the ancestor population subtree si in S, 
	replacing them by a single leaf labelled Pj, and replace labels of any leaf 
	of G belonging to species pjk with a label Pj.
	# This option minimizes the artefactual multiplication of leaves assigned to one population
	# (i.e. gene copies) when not included in collpased bush clades.
		      
	- 'collapseALLinSandG': 
	same as above, generalized to ALL populations in S : 
	collapse all Pi in S regardless of the need to collapse bush clades in G.
	Bush clades or single leaves in G that belong to 
	those collpased species populations are replaced
	by the corresponding population label.
	# This option minimizes the effective size of S tree for later reconciliation inference
	# but lead to maximal impact of artefactual multiplication of same-population leaves
	          
	- 'collapseCCinG-collapsePOPinSnotinGbutinCC': 
	same as above, generalized to ALL populations in S 
	but those present in G not as a bush ancestor : 
	collapse all Pi in S that are not represented in G, 
	and also collapse those that represent an ancestor 
	to a collpased bush clade in G.
	Bush clades or single leaves in G that belong to 
	those collpased species populations are replaced
	by the corresponding population label.
	# This option is a compromize between minimizing the effective size of S tree for later
	# reconciliation and the artefactual multiplication of gene copies
	"""
	
	def replaceCCwithLabel(cla, speorpop, fam, dold2newname, tag='CC'):
		"""collapsed clade is replaced by a single leaf with a species or species population label
		
		new G leaf label has the following structure: POPNAME_FAM_XX1234, 
		with XX=CC referring to a collapsed 'bush' clade
		or XX=SG referring to a single gene sequence (out of a collapsed 'bush' clade).
		"""
		newname = "%s_%s_%s%03d"%(speorpop, fam, tag, int(cla.split('clade')[1]))
		dold2newname[cla] = newname
		if verbose: print 'rename leaf:', newname
		
	
	def replaceCCwithSubtree(cla, ancpop, fam, dold2newname, spetree, tag='RC', restrictPopToLeaves=[]):
		"""collapsed clade is replaced by a subtree"""
		popsubtree = spetree[ancpop].deepcopybelow()
		if restrictPopToLeaves:
			popsubtree.restrictToLeaves(restrictPopToLeaves)
		if dccmaxlen:
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
		for popleaf in popsubtree:
			# new G leaf labels have the following structure: SPENAME_FAM_RC1234, with CC referring to a replaced clade
			newspename = "%s_%s_%s%03d"%(popleaf.label(), fam, tag, int(cla.split('clade')[1]))
			popleaf.edit_label(newspename)
		if verbose: print 'replace with subtree:', popsubtree
		dold2newname[cla] = tree2toBioPhylo(popsubtree) # already formats the tree in Bio.Phylo format to later integration to the parsed gene tree chain	
	
	
	if verbose: print 'mapPop2GeneTree()'
	methods = ['collapseCCinG', 'replaceCCinGasinS', 'replaceCCinGasinS-collapsePOPinSnotinG', \
	           'collapseCCinG-collapsematchingPOPinS', 'collapseALLinSandG', 'collapseCCinG-collapsePOPinSnotinGbutinCC']
	assert method in methods
	dnamepops=dict(lnamepops)
	# get file pre/sufix
	dirgt = os.path.dirname(nfingtchain1)
	bngt, extgt = os.path.basename(nfingtchain1).split('.', 1)
	#~ outbn = bngt.replace("-collapsed", '-'+method)
	outbn = bngt
	fam = bngt.rsplit('-', 1)[0]
	
	nfoutcolGtrees = "%s/%s-Gtrees.nwk"%(dirout, outbn)
	if ((os.path.exists(nfoutcolGtrees)) and (reuseOutput==2)):
		# assume work was done before, skip
		print "# (reused all) %s"%nfoutcolGtrees
		return None
	nfoutcolStree = "%s/%s-Stree.nwk"%(dirout, outbn)
	nfphyloprof = "%s/%s/%s.%s"%(dircons, phyloproftag, fam, phyloproftag)
	nfoutreflab = "%s/%s-leaflabels_Spe2Pop.txt"%(dirout, outbn)
	nfoutlabpoplab = "%s/%s-replacedlabelsbyPop.txt"%(dirout, outbn)
	nfoutfreqpopdistr = "%s/%s-PopFreqDistrib.txt"%(dirout, outbn)
	nfccmaxlen = "%s/%s/%s.%s"%(dircons, ccmaxlentag, fam, ccmaxlentag)
	
	## define mapping from the leaf set
	# try and find the consensus gene tree
	nfcongt = '.'.join([os.path.join(dirgt, bngt), extgt.replace(chain1ext, contreext)])
	nfali = '.'.join([os.path.join(dircons, collapsalntag, bngt), extgt.replace(chain1ext, aliext)])
	if os.path.exists(nfcongt):
		congt = tree2.read_nexus(nfcongt, returnDict=False, allLower=False)[0]
		gtleaflabels = congt.get_leaf_labels()
	elif os.path.exists(nfali):
		ali = AlignIO.read(nfali, 'nexus')
		gtleaflabels = [aliseq.id for aliseq in ali]
	else:
		raise IOError, "could not access either of consensus tree '%s' or alignment '%s' files"%(nfcongt, nfali)
	sspeG = set(getSpeFromGenes(gtleaflabels, species_sep=species_sep))
	genetreepopset = set([])
	for spe in sspeG:
		popname = dspe2pop.get(spe) # return None when original leaf label is not a canonical species, i.e. is a collapsed clade labelled 'cladeX'
		if popname:
			pop = dnamepops[popname]
			if len(pop)>1: genetreepopset.add(popname) # only consider non-trivial populations, i.e. grouping >= 2 species
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
	if constraintclades and os.path.exists(nfccmaxlen):
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
			dccmaxlen[cla] = ccsullsubtree.max_leaf_distance()
		try:
			with open(nfccmaxlen, 'w') as fccmaxlen:
				fccmaxlen.write('\n'.join(['%s\t%f'%cml for cml in dccmaxlen.iteritems()])+'\n')
		except IOError, e:
			sys.stderr.write("failed attempt to write max subtree lengths of the constrained clades to file '%s'\n"%nfccmaxlen)
		
	# record collapsed bush clades' leaf labels to edit in G
	colapsedcladespopset = set([])
	lancpopnames = []
	dold2newname = {}
	speciestoprune = set([])
	foutfreqpopdistr = open(nfoutfreqpopdistr, 'w') if ((not os.path.exists(nfoutfreqpopdistr)) or (not reuseOutput)) else None
	for cla, leaves in constraintclades.iteritems():
		if verbose: print '##', cla
		lspeCC = getSpeFromGenes(leaves, species_sep=species_sep)
		# identify which species population is considered the ancestor of the collapsed clades
		ancpop, ancnbcopies, popdistr = getCladeAncestor(lspeCC, spetree, dspe2pop, dnamepops, method='max_meanfreqpop', dpopnd=dpopnd, verbose=verbose)
		if foutfreqpopdistr: foutfreqpopdistr.write('\t'.join([cla]+['%s:%.3f'%(popname, popfreq) for popname, popfreq, indcount in popdistr])+'\n')
		# record the population considered ancestral to the collasped clade,
		# i.e. the one used to replace the collasped clade in G (with a pop-labeled leaf or a S-tree-copied subtree) ;
		colapsedcladespopset.add(ancpop)
		lancpopnames.append((ancpop, leaves))
		# for 'collapseCCinG*- methods, this population will appear collapsed in S;
		# for the 'replaceCCinG*' method, that will depend on representation of member species in the gene tree
		# if there are more than one, those will be chosen as representative of the population in a prune version of the Species tree
		# species to prune off the family-specific species tree are stored
			
		if method.split('-')[0].startswith('replace'):
			# evaluate presence of member species out of the clades to be colapsed
			outccancpopmembers = sspeG & set(dnamepops[ancpop])
			if len(outccancpopmembers)<=1:
				# if there are 0 or 1 members of the CC ancestral populations represented elsewhere in the Gene tree
				# better represent the CC with a single population label (CC tag)
				#if there is 1, outlier sequences of this member species will be equated with the population (SG tags, done below)
				replaceCCwithLabel(cla, ancpop, fam, dold2newname)
			else:
				# collapsed clade is replaced by a subtree covering the member species present out of the CC
				replaceCCwithSubtree(cla, ancpop, fam, dold2newname, spetree, tag='RC', restrictPopToLeaves=outccancpopmembers)
				# this population is not anymore to be collapsed in the Species tree
				colapsedcladespopset.discard(ancpop)
				speciestoprune |= (set(dnamepops[ancpop]) - outccancpopmembers)
			
		#~ elif method.split('-')[0].startswith('collapse'):
		else:
			# collapsed clade is replaced by a single leaf
			replaceCCwithLabel(cla, ancpop, fam, dold2newname, tag='CC')
			
	if foutfreqpopdistr: foutfreqpopdistr.close()
	if verbose: print '#', fam, 'colapsedcladespopset', colapsedcladespopset
	# record the set of populations absent in G;
	# will appear collapsed in S in methods '*-collapsePOPinSnotinG'
	notingttreepopset = (set(dnamepops.keys()) - genetreepopset)
	# record the set of populations absent in G unless they are considered ancestral to any collasped clades;
	# will appear collapsed in S in methods '*-collapsePOPinSnotinGbutinCC'
	notingtbutincctreepopset = notingttreepopset | colapsedcladespopset
	
	if verbose: print '#', fam, 'dold2newname (CC)', dold2newname
	## prepare name reference dict for editing labels in gene tree sample
	dpop2SGnggen = {popname:0 for popname in dnamepops}
	if 'inS' in method:
		# collapse populations in S
		if method=='collapseALLinSandG':
			lnamepops2collapse = lnamepops
		elif method.endswith('collapsePOPinSnotinG'):
			lnamepops2collapse = [(name, pop) for name, pop in lnamepops if (name in notingttreepopset)]
		elif method=='collapseCCinG-collapsematchingPOPinS':
			lnamepops2collapse = [(name, pop) for name, pop in lnamepops if (name in colapsedcladespopset)]
		elif method.endswith('collapsePOPinSnotinGbutinCC'):
			lnamepops2collapse = [(name, pop) for name, pop in lnamepops if (name in notingtbutincctreepopset)]
		# S tree with some populations collapsed as ancestors and some not, used later for reconciliation inference on this particular family
		if verbose: print '#', fam, 'lnamepops2collapse:\n %s'%repr_long_list(lnamepops2collapse)
		dspe2pop2collapse = getdspe2pop(lnamepops2collapse)
		try:
			getsavepoptree(nfoutcolStree, spetree=spetree, lnamepops2collapse=lnamepops2collapse, speciestoprune=speciestoprune, verbose=max(0, verbose-1))
		except ValueError, e:
			sys.stderr.write(str(e)+'\n')
			sys.stderr.write("Error processing family %s\n"%fam)
			return 1
		# then record other single leaf labels belonging to collapsed pops in S to edit in G
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
				if isinstance(newlaborst, str):
					newlabs = [newlaborst]
				elif isinstance(newlaborst, tree2.Node):
					newlabs = newlaborst.get_leaf_labels()
				elif isinstance(newlaborst, BaseTree.TreeElement):
					newlabs = [tip.name for tip in newlaborst.get_terminals()]
				for newlab in newlabs:
					newlabpop, fam, ngtag = newlab.split('_')
					dpoptyperepcount.setdefault(newlabpop, {rt:0 for rt in reptypes})[ngtag[:2]] += 1
					foutreflab.write('\t'.join((oldlab, newlab))+'\n')
	if dpoptyperepcount or ((not os.path.exists(nfoutlabpoplab)) or (not reuseOutput)):
		with open(nfoutlabpoplab, 'w') as foutlabpoplab:
			foutlabpoplab.write('\t'.join(['Pop_or_Spe_id']+reptypes)+'\n')
			for pop in dpoptyperepcount:
				foutlabpoplab.write('\t'.join([pop]+[str(dpoptyperepcount[pop][reptype]) for reptype in reptypes])+'\n')
			
	## apply to gene tree sample
	# load gene tree sample from Newick concatenate tree file
	if ((not os.path.exists(nfoutcolGtrees)) or (not reuseOutput)):
		nfingtchains = [nfingtchain1.replace(chain1ext, chain1ext.replace('1', str(k))) for k in range(1, nbchains+1)]
		parseChain(nfingtchains, dold2newname, nfchainout=nfoutcolGtrees, verbose=verbose)
	else:
		print "# (reused) %s"%nfoutcolGtrees
		
def usage():
	s =  'For population assignation:\n'
	s += ' python replace_species_by_pop_in_gene_trees.py -S /path/to/species_tree\n'
	s += 'For replacement of species labels with population labels in gene tree samples:\n'
	s += ' python replace_species_by_pop_in_gene_trees.py -S /path/to/[TIMED.]species.tree -G /path/to/list.of.gene.tree.sample.chain1.files -c /path/to/folder.of.collapsed.clade.info.files -o /path/to/output.folder \\\n'
	s += ' [--populations=/path/to/species2population_map --population_tree=/path/to/population_annotated_species_tree] [OTHER OPTIONS]\n'
	s += ' Options:\n'
	s += '  Input file options:\n'
	s += '  --populations\t\t\tpath to file with already computed mapping from species populations labels to species leaf labels.\n'
	s += '  --population_tree\t\t\tpath to Newick-formatted tree file of an already computed mapping of populations ancestors onto the species tree.\n'
	s += '  --population_node_distance\t\tpath to pre-computed matrix of inter-node distances between populations in the species tree.\n'
	s += '  --dir_full_gene_trees\t\tpath to folder of original (before clade collapsing) gene trees (can be used to scale subtree grafts)\n\t\t\t\t(overridden by the presence of a adequate file in the collapsed.clade.info.files subfolder \'max_subtree_lengths\').\n'
	s += '  --pop_stem_conds\t\tconditions on clade stem to select populations in the species tree (see mark_unresolved_clades.py)\n'
	s += '  --within_pop_conds\t\tconditions on clade subtree to select populations in the species tree (see mark_unresolved_clades.py)\n'
	s += '  --'+mapPop2GeneTree.__doc__+'\n'
	s += '  --threads\t\t\tnumber of parralel processes to run.\n'
	s += '  --reuse\t\t[0,1,2]\tif 0: normal behaviour;\n'
	s += '	       \t\t\t\tif 1: skip computing and writing output files which name matches already existing file;\n'
	s += '	       \t\t\t\t      BEWARE! This can save a lot of time but may produce inconsistent output file sets!\n'
	s += '	       \t\t\t\tif 2: skip processing the entire gene tree chain if the expected output gene tree chain file (with labels replaced) exists;\n'
	s += '  --verbose={0,1,2}\tverbose mode, from none to plenty.\n'
	s += '  -v\t\t\tequivalent to --verbose=1.\n'
	return s

	
if __name__=='__main__':	
	opts, args = getopt.getopt(sys.argv[1:], 'G:c:S:o:hv', ['dir_out=', 'method=', \
															'populations=', 'population_tree=', 'population_node_distance=', 'dir_full_gene_trees=', \
															'pop_stem_conds=', 'within_pop_conds=', \
															'threads=', 'reuse=', 'verbose=', 'help'])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	nfspetree = dopt['-S']
	nflnfingtchain = dopt.get('-G')
	dircons = dopt.get('-c')
	dirout = dopt.get('-o')
	methods = dopt.get('--method', 'collapseCCinG,replaceCCinGasinS,replaceCCinGasinS-collapsePOPinSnotinG').split(',')
	nfpops = dopt.get('--populations')
	nfpoptree = dopt.get('--population_tree')
	nfdpopnd = dopt.get('--population_node_distance')
	dirfullgt = dopt.get('--dir_full_gene_trees')
	psc = eval(dopt.get('--pop_stem_conds', 'None'))
	wpc = eval(dopt.get('--within_pop_conds', 'None'))
	nbthreads = int(dopt.get('--threads', 1))
	reuseOutput = int(dopt.get('--reuse', 0))
	verbose = int(dopt.get('--verbose', ('-v' in dopt)))
	
	### define species populations
	spetree, dspe2pop, lnamepops, dpopnd = inferPopfromSpeTree(nfspetree, populations=nfpops, \
                                                               nfpopulationtree=nfpoptree, nfdpopnodedist=nfdpopnd, \
                                                               pop_stem_conds=psc, within_pop_conds=wpc, \
                                                               nbthreads=nbthreads, verbose=verbose)
				
	if nflnfingtchain:
	
		if not (dircons and dirout): raise ValueError, "must provide value for options -c and -o when passing gene tree samples to -G"
		ldircreate = [os.path.join(dircons, phyloproftag)]+[os.path.join(dirout, meth) for meth in methods]
		for d in ldircreate:
			if not os.path.isdir(d):
				os.mkdir(d)
		
		### map populations on gene trees	
		with open(nflnfingtchain, 'r') as flnfingtchain:
			lnfingtchain = [line.strip('\n') for line in flnfingtchain]
		
		def mapPop2GeneTreeSetArgs(nfingtchain):
			print '#', nfingtchain
			for method in methods:
				print 'method:', method
				try:
					mapPop2GeneTree(nfingtchain, dircons, os.path.join(dirout, method), method, spetree, dspe2pop, lnamepops, dpopnd, \
									reuseOutput=reuseOutput, phyloproftag=phyloproftag, dirfullgt=dirfullgt, verbose=verbose)
				except Exception, e:
					if type(e) is RuntimeError:
						print "met RuntimeError on file '%s':"%nfingtchain, e
						if str(e).startswith('maximum recursion depth exceeded'):
							if sys.getrecursionlimit()<10000:
								sys.setrecursionlimit(sys.getrecursionlimit()+2000)
								print "re-run with increased recursion depth +2000"
								mapPop2GeneTreeSetArgs(nfingtchain)
							else:
								print "Error: failed to run on file '%s':"%nfingtchain
								return 1
					else:	
						print 'Caught exception:'
						traceback.print_exc()
						raise e
			print '# # # #\n'
		
		if nbthreads==1:
			for nfingtchain in lnfingtchain:
				if verbose: print nfingtchain
				mapPop2GeneTreeSetArgs(nfingtchain)
		else:
			pool = mp.Pool(processes=nbthreads)
			res = pool.map(mapPop2GeneTreeSetArgs, lnfingtchain)
			#~ res = [pool.apply_async(mapPop2GeneTreeSetArgs, nfingtchain) for nfingtchain in lnfingtchain] #, chunksize=5

