#!/usr/bin/python

import tree2
import sys, os, getopt, glob
from mark_unresolved_clades import select_clades_on_conditions, median
import copy

import multiprocessing as mp
from Bio import File, AlignIO, Align
from Bio.Phylo import BaseTree, NewickIO, NexusIO
from Bio.Phylo import _io as PhyloIO

supported_formats = {'newick': NewickIO, 'nexus': NexusIO}

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

def parseChain(lnfchains, dold2newname, nfchainout=None, inchainfmt='nexus', outchainfmt='newick'):
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
			try:
				tree = dhandles[k].next()
			except StopIteration:
				sometreeleft = False
				# stops for all chains here, even if they had different length (whichc shoulld not be)
				break # the for k loop
			if dold2newname:
				for tip in tree.get_terminals():
					nutipname = dold2newname.get(tip.name)
					if nutipname: tip.name = nutipname
			treebuffer.append(tree)
			ntree += 1
			if len(treebuffer) >= buffsize:
				if ntree <= buffsize:
					PhyloIO.write(treebuffer, nfout, outchainfmt)
				else:
					treeappend(treebuffer, nfout, outchainfmt)
				treebuffer = []
	if treebuffer: treeappend(treebuffer, nfout, outchainfmt)
	print nfout, '...done (%d trees)'%ntree
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
	lpops = select_clades_on_conditions(spetree , clade_stem_conds=pop_stem_conds, within_clade_conds=within_pop_conds, nested=nested, inclusive=inclusive, **kw)
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

def collapsePopulationInSpeciesTree(spetree, lnamepops, nested=True, keepUltrametric=True, verbose=0):
	"""if population are nested, the order of the input population (name, (members,)) list matters!
	verbose=0: silent, verbose=1: general progression, verbose=2: detail of tree edition.
	"""
	# first mark the population ancestor node
	spetree.complete_node_ids(order=2) # post-order traversal numeration
	poptree = copy.deepcopy(spetree)
	dnamepop = dict(lnamepops)
	# list of nodes independent of tree structure (i.e. pruning or unlinking nodes in the tree won't affect list sequence]
	lnamepopancs = [(popname, poptree.mrca(popspe)) for popname, popspe in lnamepops]
	# ensure sorting pop. ancestor nodes following post-order traversal numeration in the original tree
	lnamepopancs.sort(key=lambda x: x[1].nodeid())
	lpopname, lpopancs = zip(*lnamepopancs)
	lpopancids = [popanc.nodeid() for popanc in lpopancs] # increasing sequence of node ids
	for a, (popname, popanc) in enumerate(lnamepopancs):
		lpopspe = dnamepop[popname]
		if verbose: print popname
		if len(lpopspe)==1: 
			# leaf/single-species population; skip
			continue
		sublen = popanc.max_leaf_distance()
		if nested:
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
							# mark the population name on that leaf; allow for extra branch langth to be kept relative to the true population ancestor 
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
				if keepUltrametric and sublen!='NA': popanc += sublen	# reprots ditances to longest leaf lineage below the ancestor so that an ultrametric tree keeps being ultrametric
				if verbose: print 'as_leaf', popname, '+=', sublen
						
		else:
			# assume populations are monophyletic and the clade can be removed at once
			assert set(popanc.get_leaf_labels())==set(lpopspe)
			popanc.as_leaf(newlabel=popname, silent=1-(verbose-1))
	return poptree

#~ nfspetree = '/home/flassall/PanteroDB/core_genome/raxml_tree/RAxML_bipartitions.concat_core_cds_464entero.rooted.codes.nwk'
#~ spetree = tree2.AnnotatedNode(file=nfspetree)
#~ spetree.complete_internal_labels(order=0, ffel=True) # pre-order traversal labeling
#~ lnamepops = speciesTreePopulations(spetree, pop_stem_conds=[('lg', '>=', 0.0002), ('bs', '>=', 80)], within_pop_conds=[('max', 'lg', '<=', 0.0002, -1)], testLeaves=True, nested=True, inclusive=False, taglen=3)
#~ dspe2pop = getdspe2pop(lnamepops, spetree)
#~ poptree = collapsePopulationInSpeciesTree(spetree, lnamepops, nested=True)
#~ poptree.figtree()
#~ annotatePopulationInSpeciesTree(spetree, lnamepops)

def countPops(dspecount, spetree, dspe2pop):
	dpopcount = {}
	for spe in dspecount:
		pop = dspe2pop[spe]	# have to make sure before all species are included in a population
		dpopcount[pop] = dpopcount.setdefault(pop, 0) + dspecount[spe]
	popcount = sorted(dpopcount.items(), key=lambda x: x[1])
	return popcount

def meanFreqPops(dspecount, spetree, dspe2pop, dnamepops):
	dpopcount = {}
	for spe in dspecount:
		popname = dspe2pop[spe]	# have to make sure before all species are included in a population
		dpopcount[popname] = dpopcount.setdefault(popname, 0) + dspecount[spe]
	meanpopfreq = [(popname, float(popcount)/len(dnamepops[popname])) for popname, popcount in dpopcount.items()]
	return sorted(meanpopfreq, key=lambda x: x[1])

def medianFreqPops(dspecount, spetree, dspe2pop, dnamepops):
	dpopcount = {}
	for spe in dspecount:
		popname = dspe2pop[spe]	# make sure all species are included in a population
		dpopcount.setdefault(popname, []).append(dspecount[spe])
	medpopfreq = [(popname, median(popcount)) for popname, popcount in dpopcount.items()]
	return sorted(medpopfreq, key=lambda x: x[1])

def getSpeFromGenes(leaflabs, species_sep='_', **kw):
	return [leaflab.split(species_sep)[0] for leaflab in leaflabs]

def getCladeAncestor(lspe, spetree, method='max_meanfreqpop', **kw):
	## first a census of represented species
	dspecount = {spe:lspe.count(spe) for spe in set(lspe)}
	## then use that to determine who was there first
	crit, quant = method.split('_')
	# rank populations according to targeted quantity
	if quant=='countpop':
		popdistr = countPops(dspecount, spetree, **kw)
	elif quant=='meanfreqpop':
		popdistr = meanFreqPops(dspecount, spetree, **kw)
	elif quant=='medianfreqpop':
		popdistr = medianFreqPops(dspecount, spetree, **kw)
	# and select given optimal criterion
	if crit=='max':
		ancpop, popquant = popdistr[-1]
	elif crit=='min':
		ancpop, popquant = popdistr[0]
	if quant.endswith('freqpop'):
		# keep track of roughly how many gene copies were there in the clade
		# to avoid biasing to much the reconciliation if there was 
		# a conserved ancestral duplication above the clade's node
		ancnbcopies = popquant
	else:
		ancnbcopies = 1
	return (ancpop, ancnbcopies, popdistr)
	

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

def getsavepoptree(nfpoptreeout, poptree=None, spetree=None, lnamepops2collapse=None, verbose=False):
	if not poptree: poptree = collapsePopulationInSpeciesTree(spetree, lnamepops2collapse, nested=True, verbose=verbose)
	# save along the output gene tree sample file (species pruning being family-specific)
	cpoptree = copy.deepcopy(poptree)
	for node in cpoptree: # for compatibility with ALE
		if not node.is_leaf(): node.edit_label('')
	cpoptree.write_newick(nfpoptreeout)
	return poptree

def inferPopfromSpeTree(nfspetree, \
         populations=None, populationtree=None, \
         pop_stem_conds=None, within_pop_conds=None, \
         verbose=False):
	
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
		if populationtree:
			poptree = tree2.AnnotatedNode(file=populationtree)
			if verbose: print 'load population tree from \'%s\''% populationtree
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
				
	poptree = getsavepoptree("%s_collapsedPopulations.nwk"%(bnst), poptree=poptree, spetree=spetree, lnamepops2collapse=lnamepops, verbose=verbose)
	
	## annotate ancestor nodes of populations on S
	annotatePopulationInSpeciesTree(spetree, lnamepops)
	dspe2pop = getdspe2pop(lnamepops, spetree)
	
	return (spetree, dspe2pop, lnamepops)
	
def mapPop2GeneTree(nfingtchain1, dircons, spetree, dspe2pop, lnamepops, dirout, method, reuseOutput=0, \
                    chain1ext='mb.run1.t', nbchains=2, aliext='codes.nex', contreext='mb.con.tre', constrainttag='mbconstraints', \
                    species_sep='_', verbose=False, **kw):
	"""method can be 'replaceGasinS', 'collapseinSandG' or 'collapseALLinSandG'
	
	P is the set of populations {Pi}, each being a set of member species leaves (pik) from the species tree S, mappping to a subtree si.
	C is the set of "bush" clades {Ci} in the gene tree G, defined by a well-suported stem and no strong support for its topology within
	with Pi the population identified as the ancestor of a "bush" clade Ci,
	- 'collapseinG': will keep S as is, and replace Ci in G by a single leaf labelled Pi
	
	(!!! the following methods lead to change of G-to-S mapping cardinality, i.e. some leaves representing different species assigned to a same population 
	     will artefactually be seen as multiple gene copies for one single species population)
	- 'collapseinSandG': collapse both the population subtrees si in S and Ci in G, replacing them by a single leaf labelled Pi, 
						 and replace label of any leaf of G belonging to species pik with a label Pi
		      # This option minimizes the artefactual multiplication of leaves assigned to one population (i.e. gene copies) when not included in collpased bush clades.
		      
	- 'collapseALLinSandG': same as above, generalized to ALL populations in S : 
	                        collapse all Pi in S regardless of the need to collapse bush clades in G, and replace matching bushes/leaves in G 
	          # This option minimizes the effective size of S tree for later reconciliation inference
	          
	- 'collapseALLinSbutinG': same as above, generalized to ALL populations in S but those present in G not as a bush ancestor : 
	                             collapse all Pi in S that not represented in G or represent a collpased bush clade in G, and replace matching bushes/leaves in G 
	          # This option is a compromize between minimizing the effective size of S tree for later reconciliation and the artefactual multiplication of gene copies
	"""
	if verbose: print 'mapPop2GeneTree()'
	methods = ['collapseinG', 'collapseinSandG', 'collapseALLinSandG', 'collapseALLinSbutinG']
	assert method in methods
	dnamepops=dict(lnamepops)
	# get file pre/sufix
	dirgt = os.path.dirname(nfingtchain1)
	bngt, extgt = os.path.basename(nfingtchain1).split('.', 1)
	outbn = bngt.replace("-collapsed", '-'+method)
	fam = bngt.rsplit('-', 1)[0]
	
	nfoutcolGtrees = "%s/%s-Gtrees.nwk"%(dirout, outbn)
	if ((os.path.exists(nfoutcolGtrees)) and (reuseOutput==2)):
		# assume work was done before, skip
		print "# (reused all) %s"%nfoutcolGtrees
		return None
	nfoutcolStree = "%s/%s-Stree.nwk"%(dirout, outbn)
	nfphyloprof = "%s/%s.phyloprofile"%(dircons, fam)
	nfoutreflab = "%s/%s-leaflabels_Spe2Pop.txt"%(dirout, outbn)
	nfoutlabpoplab = "%s/%s-replacedlabelsbyPop.txt"%(dirout, outbn)
	nfoutfreqpopdistr = "%s/%s-PopFreqDistrib.txt"%(dirout, outbn)
	
	## define mapping from the leaf set
	# try and find the consensus gene tree
	nfcongt = '.'.join([os.path.join(dirgt, bngt), extgt.replace(chain1ext, contreext)])
	nfali = '.'.join([os.path.join(dircons, bngt), extgt.replace(chain1ext, aliext)])
	if os.path.exists(nfcongt):
		congt = tree2.read_nexus(nfcongt, returnDict=False, allLower=False)[0]
		gtleaflabels = congt.get_leaf_labels()
	elif os.path.exists(nfali):
		ali = AlignIO.read(nfali, 'nexus')
		gtleaflabels = [aliseq.id for aliseq in ali]
	else:
		raise IOError, "could not access either of consensus tree '%s' or alignment '%s' files"%(nfcongt, nfali)
	treepopset = set([])
	for spe in getSpeFromGenes(gtleaflabels, species_sep=species_sep):
		popname = dspe2pop.get(spe)
		pop = dnamepops.get(popname, []) # return empty pop when original leaf label is not a canonical species, i.e. is a collapsed clade labelled 'cladeX'
		if len(pop)>1: treepopset.add(popname) # only consider non-trivial populations, i.e. grouping >= 2 species
	# load list of gene subtree constrained clades = clades collapsed in the full gene tree G
	if verbose: print "get collapsed clades from files matching '%s/%s*%s*' pattern"%(dircons, fam, constrainttag)
	lnfcons = glob.glob("%s/%s*%s*"%(dircons, fam, constrainttag))
	if verbose: print lnfcons
	constraintclades = parseMrBayesConstraints(lnfcons)
	# write out the mapping of clades to species as phyloprofiles based on the species tree
	if ((not os.path.exists(nfphyloprof)) or (not reuseOutput)):
		famprof = get_clade_phyloprofiles(constraintclades, spetree)
		write_out_clade_phyloprofiles({fam:famprof}, nfphyloprof)
	# record collapsed bush clades' leaf labels to edit in G
	colpasedcladespopset = set([])
	lancpopnames = []
	dold2newname = {}
	foutfreqpopdistr = open(nfoutfreqpopdistr, 'w') if ((not os.path.exists(nfoutfreqpopdistr)) or (not reuseOutput)) else None
	for cla, leaves in constraintclades.iteritems():
		lspe = getSpeFromGenes(leaves, species_sep=species_sep)
		#~ clapopset = set([dspe2pop[spe] for spe in lspe])
		#~ colpasedcladespopset |= clapopset
		# identify which species population is considered the ancestor of the collapsed clades
		ancpop, ancnbcopies, popdistr = getCladeAncestor(lspe, spetree, method='max_meanfreqpop', dspe2pop=dspe2pop, dnamepops=dnamepops)
		# new G leaf label has the following structure: POPNAME_C1234, with CC referring to a collapsed clade
		newname = "%s_%s_CC%03d"%(ancpop, fam, int(cla.split('clade')[1]))
		colpasedcladespopset.add(ancpop)
		dold2newname[cla] = newname
		lancpopnames.append((ancpop, leaves))
		if foutfreqpopdistr: foutfreqpopdistr.write('\t'.join([cla]+['%s:%.3f'%(pop, freq) for pop, freq in popdistr])+'\n')
	if foutfreqpopdistr: foutfreqpopdistr.close()
	# record the set of populations considered ancestral to the collasped clades, i.e. those used to replace the label of the collasped clades in G ; will appear collapsed in S for all methods
	if verbose: print '#', fam, 'colpasedcladespopset', colpasedcladespopset
	# record the set of populations absent in G unless they are considered ancestral to any collasped clades; will appear collapsed in S in method 'collapseALLinSbutinandG'
	notingtbutincctreepopset = (set(dnamepops.keys()) - treepopset) | colpasedcladespopset
	
	if verbose: print '#', fam, 'dold2newname (CC)', dold2newname
	## prepare name reference dict for editing labels in gene tree sample
	dpop2SGnggen = {popname:0 for popname in dnamepops}
	if 'inS' in method:
		# collapse populations in S
		if method=='collapseinSandG':
			lnamepops2collapse = [(name, pop) for name, pop in lnamepops if (name in colpasedcladespopset)]
		elif method=='collapseALLinSandG':
			lnamepops2collapse = lnamepops
		elif method=='collapseALLinSbutinG':
			lnamepops2collapse = [(name, pop) for name, pop in lnamepops if (name in notingtbutincctreepopset)]
		# S tree with some populations collapsed as ancestors and some not, used later for reconciliation inference on this particular family
		if verbose: print '#', fam, 'lnamepops2collapse', lnamepops2collapse
		dspe2pop2collapse = getdspe2pop(lnamepops2collapse)
		try:
			getsavepoptree(nfoutcolStree, spetree=spetree, lnamepops2collapse=lnamepops2collapse, verbose=verbose)
		except ValueError, e:
			sys.stderr.write(str(e)+'\n')
			sys.stderr.write("Error processing family %s\n"%fam)
			return 1
		# then record other leaf labels belonging to collapsed pops in S to edit in G
		for leaflab in gtleaflabels:
			spelab = leaflab.split('_', 1)[0]
			poplab = dspe2pop2collapse.get(spelab)
			if poplab and poplab!=spelab:
				# new G leaf label has the following structure: POPNAME_C1234, with SG referring to a single gene sequence
				newname = "%s_%s_SG%03d"%(poplab, fam, dpop2SGnggen[poplab])
				dpop2SGnggen[poplab] += 1
				dold2newname[leaflab] = newname
				
	if verbose: print '#', fam, 'dold2newname (CC+SG)', dold2newname
	## output label transaltion reference
	dpoptyperepcount = {}
	if ((not os.path.exists(nfoutreflab)) or (not reuseOutput)):
		with open(nfoutreflab, 'w') as foutreflab:
			for oldlab, newlab in dold2newname.iteritems():
				newlabpop, fam, ngtag = newlab.split('_')
				dpoptyperepcount.setdefault(newlabpop, {'CC':0, 'SG':0})[ngtag[:2]] += 1
				foutreflab.write('\t'.join((oldlab, newlab))+'\n')
	if dpoptyperepcount or ((not os.path.exists(nfoutlabpoplab)) or (not reuseOutput)):
		with open(nfoutlabpoplab, 'w') as foutlabpoplab:
			foutlabpoplab.write('Pop_id\tCC\tSG\n')
			for pop in dpoptyperepcount:
				foutlabpoplab.write('%s\t%d\t%d\n'%((pop, dpoptyperepcount[pop]['CC'], dpoptyperepcount[pop]['SG'])))
			
	## apply to gene tree sample
	# load gene tree sample from Newick concatenate tree file
	if ((not os.path.exists(nfoutcolGtrees)) or (not reuseOutput)):
		nfingtchains = [nfingtchain1.replace(chain1ext, chain1ext.replace('1', str(k))) for k in range(2, nbchains)]
		parseChain(nfingtchain1, dold2newname, nfchainout=nfoutcolGtrees, sisterchains=nfingtchains)
	else:
		print "# (reused) %s"%nfoutcolGtrees
		
def usage():
	s =  'For population assignation:\n'
	s += ' python replace_species_by_pop_in_gene_trees.py -S /path/to/species_tree\n'
	s += 'For replacement of species labels with population labels in gene tree samples:\n'
	s += ' python replace_species_by_pop_in_gene_trees.py -S /path/to/[TIMED.]species.tree -G /path/to/list.of.gene.tree.sample.chain1.files -c /path/to/folder.of.collapsed.clade.info.files -o /path/to/output.folder \\\n'
	s += ' [--populations=/path/to/species2population_map --population_tree=/path/to/population_annotated_species_tree] [OTHER OPTIONS]\n'
	s += '  Options:\n'
	s += '  --populations\t\t\tpath to file with already computed mapping from species populations labels to species leaf labels.\n'
	s += '  --population_tree\t\t\tpath to Newick-formatted tree file of an already computed mapping of populations ancestors onto the species tree.\n'
	s += '  --pop_stem_conds\t\tconditions on clade stem to select populations in the species tree (see mark_unresolved_clades.py)\n'
	s += '  --within_pop_conds\t\tconditions on clade subtree to select populations in the species tree (see mark_unresolved_clades.py)\n'
	s += '  --'+mapPop2GeneTree.__doc__+'\n'
	s += '  --threads\t\t\tnumber of parralel processes to run.\n'
	s += '  --reuse\t\t[0,1,2]\tif 0: normal behaviour;\n'
	s += '	       \t\t\t\tif 1: skip computing and writing output files which name matches already existing file;\n'
	s += '	       \t\t\t\t      BEWARE! This can save a lot of time but may produce inconsistent output file sets!\n'
	s += '	       \t\t\t\tif 2: skip processing the entire gene tree chain if the expected output gene tree chain file (with labels replaced) exists;\n'
	s += '  --verbose\t\t\tverbose mode.\n'
	return s

	
if __name__=='__main__':	
	opts, args = getopt.getopt(sys.argv[1:], 'G:c:S:o:hv', ['method=', 'populations=', 'population_tree=', \
															 'pop_stem_conds=', 'within_pop_conds=', \
															 'threads=', 'reuse=', 'verbose', 'help'])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	nfspetree = dopt['-S']
	nflnfingtchain = dopt.get('-G')
	dircons = dopt.get('-c')
	dirout = dopt.get('-o')
	methods = dopt.get('--method', 'collapseinSandG,collapseALLinSandG,collapseALLinSbutinG').split(',')
	pops = dopt.get('--populations')
	poptree = dopt.get('--population_tree')
	psc = eval(dopt.get('--pop_stem_conds', 'None'))
	wpc = eval(dopt.get('--within_pop_conds', 'None'))
	nbthreads = int(dopt.get('--threads', 1))
	reuseOutput = int(dopt.get('--reuse', 0))
	verbose = ('-v' in dopt) or ('--verbose' in dopt)
	
	### define species populations
	spetree, dspe2pop, lnamepops = inferPopfromSpeTree(nfspetree, \
                                                       populations=pops, populationtree=poptree, \
                                                       pop_stem_conds=psc, within_pop_conds=wpc, \
                                                       verbose=verbose)
			
				
	if nflnfingtchain:
		if not (dircons and dirout): raise ValueError, "must provide value for options -c and -o when passig gene tree samples to -G"
		### map populations on gene trees	
		with open(nflnfingtchain, 'r') as flnfingtchain:
			lnfingtchain = [line.strip('\n') for line in flnfingtchain]
		
		def mapPop2GeneTreeSetArgs(nfingtchain):
			print '#', nfingtchain
			for method in methods:
				print 'method:', method
				mapPop2GeneTree(nfingtchain, dircons, spetree, dspe2pop, lnamepops, dirout, method, reuseOutput=reuseOutput, verbose=verbose)
			print '# # # #\n'
		
		if nbthreads==1:
			for nfingtchain in lnfingtchain:
				if verbose: print nfingtchain
				mapPop2GeneTreeSetArgs(nfingtchain)
		else:
			pool = mp.Pool(processes=nbthreads)
			res = pool.map(mapPop2GeneTreeSetArgs, lnfingtchain)
			#~ res = [pool.apply_async(mapPop2GeneTreeSetArgs, nfingtchain) for nfingtchain in lnfingtchain] #, chunksize=5

