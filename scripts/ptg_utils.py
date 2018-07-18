#!/usr/bin/python

import sys
from Bio import File
from Bio.Phylo import BaseTree, NewickIO, NexusIO, _io as PhyloIO
from StringIO import StringIO
from random import randint

supported_formats = {'newick': NewickIO, 'nexus': NexusIO}

### general purpose functions

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
	
def addStrtoTuple(a, b):
	"""addition for tuples, suporting appending of single str element; always return a single tuple"""
	if type(a) is str:
		if type(b) is str:
			return (a, b)
		elif type(b) is tuple:
			return (a,) + b
		else:
			raise TypeError, "unexpected type for 'b': %s"%repr(b)
	elif type(a) is tuple:
		if type(b) is str:
			return a + (b,)
		elif type(b) is tuple:
			return a + b
		else:
			raise TypeError, "unexpected type for 'b': %s"%repr(b)
	else:
		raise TypeError, "unexpected type for 'a': %s"%repr(a)

def mean(seq, ignoreNull=True):
	l = [float(k) for k in seq if ((k is not None) or (not ignoreNull))]
	if not l: return None
	return sum(l)/len(l)

def median(seq, ignoreNull=True):
	l = [k for k in seq if ((k is not None) or (not ignoreNull))]
	l.sort()
	L = len(l)
	if not l: return None
	nmed = L/2
	if L%2==1: return float(l[nmed])
	else: return (float(l[nmed-1])+float(l[nmed]))/2

def var(seq, correct=1):
	m = mean(seq)
	if m is None: return None
	n = len([k for k in seq if (k is not None)])-correct
	if n < 1:
		return 'inf'		
	else:
		d = [(float(x) - m)**2 for x in seq]
		return sum(d)/n

## Functions related to Species-to-Population matching and Gene Tree Clade Collpasing

def getSpeFromGenes(leaflabs, species_sep='_', **kw):
	return [leaflab.split(species_sep)[0] for leaflab in leaflabs]

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

def loadRecGeneTreeLabelAliases(nfgenefamlist, dircons=None, dirrepl=None, sep='\t', nbthreads=1, verbose=False):
	if not nfgenefamlist: 
		return []
	with open(nfgenefamlist, 'r') as fgenefamlist:
		header = fgenefamlist.readline().strip(' \n').split(sep) 
		# remove trailing whitespace present after gene_family_id ; account for PostgreSQL footer line count of fetched rows
		genefamlist = [dict(zip(header, line.replace(' '+sep, sep).strip(' \n').split(sep))) for line in fgenefamlist if (not line.endswith('rows)\n'))] 
		if verbose: print 'load restricted gene lineage and family table with header:', header, 'and %d lines'%len(genefamlist)
	dreplacedlab = {}
	if dircons or dirrepl:
		for genefam in genefamlist:
			fam = genefam['gene_family_id']
			if dircons:
				globcons = "%s/%s/*%s*"%(dircons, fam, constrainttag)
				if verbose: print globcons
				lnfcons = glob.glob(globcons)
				constraintclades = parseMrBayesConstraints(lnfcons)
			else:
				constraintclades = []
			if dirrepl:
				globreplaced = "%s/%s*%s"%(dirrepl, fam, replacedtag)
				if verbose: print globreplaced
				lnfreplaced = glob.glob(globreplaced)
				if len(lnfreplaced)!=1:
					sys.stderr.write("Warning! will use first file among those found matching pattern '%s'\nas reference for collapsed gene tree tip label aliases:\n%s\n"%( globreplaced, '\n'.join(lnfreplaced) ))
				nfreplaced = lnfreplaced[0]
				prevclaorgene = ''
				with open(nfreplaced, 'r') as freplaced:
					for line in freplaced:
						claorgene, repllab = line.rstrip('\n').split('\t')
						if claorgene in constraintclades:
							for genelab in constraintclades[claorgene]:
								dreplacedlab[genelab] = repllab
						else:
							if claorgene == prevclaorgene: continue
							# NOTE: this will always leave only one leaf label associated with a collapsed clade,
							# even when gene tree collapsed clade was replaced by replacement clde (tagged 'RC')
							# with multiple leaves, mimicking the species tree : because these clades are artificially
							# introduced in the gene trees and will always show perfect agreement with the species trere,
							# spurious association of S events with total support would be seen between similar RC clades
							# in different gene families.
							# (events recorded withtin such RC clade are actually discarded in pAr.parseRecGeneTree(),
							# as they detected as clade encompassing all leaves with identical RC tag)
							dreplacedlab[claorgene] = repllab
						prevclaorgene = claorgene
	
	for genefam in genefamlist:
		genelab = genefam.get('cds_code')
		if genelab in dreplacedlab:
			genefam['replaced_cds_code'] = dreplacedlab[genelab]
	return genefamlist

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

#~ def colorWheel(i, n):
	#~ assert i < n
	#~ j = float(i)
	#~ return tuple(int((256*((float(k)/3)+(j/n)))%256) for k in range(3))
def colorWheel():
	return tuple(randint(0, 255) for k in range(3))

def colour_tree_parts(tree, **kw):
	"""paints a tree with a colour pallette; take a tree and a list of leaf sets as input.
	
	the painting method is defined by the name of argument the list of leaf sets is passed with:
	'constraints' defines contraint clades:
		Paints the clades (i.e. all the branches at the stema and below a node) 
		in a reverse order relative to input leaf set list, as it is supposed 
		to be the ouput of select_clades_on_conditions(..., order=2), explored 
		from tips to root (post-order traversal), gradually removing clades, 
		which means ealier picked constrained clades can be nested in later ones. 
		The later clades are thus painted fist, and then nested ealier clades are 
		painted over.
	'leafgroups' defines arbitrary groups of leaves
	
	other optional keyword arguments:
	'palette' defines the colour palette; if not specified, defaults to a random RGB combination.
	'dalias' provides an optional alias dict to search for leaves in the tree.
	"""
	dalias = kw.get('dalias', {})
	if ('constraints' in kw) and ('leafgroups' in kw):
		raise ValueError, "cannot pass both arguments 'constraints' and 'leafgroups' at the same time"
	leafsets = kw.get('constraints', kw.get('leafgroups'))
	cols = kw.get('palette', [colorWheel() for i in range(len(leafsets))])	
	if ('constraints' in kw):
		for k, cons in enumerate(reversed(leafsets)):
			# start from the last assigned node constraint, so that ealier nested clades will be painted over at a further iteration
			n = tree.mrca([dalias.get(ll, ll) for ll in cons])
			for c in n.get_all_children():
				c.edit_color(cols[k])
	elif ('leafgroups' in kw):
		for k, lg in enumerate(leafsets):
			for ll in lg:
				tree[dalias.get(ll, ll)].edit_color(cols[k])
	else:
		raise ValueError, "need to pass on leaf set list through either 'constraints' or 'leafgroups' arguments"

def colour_tree_with_constrained_clades(tree, constraints, **kw):
	"""paints constrained clades with a colour pallette; take a tree and a list of leaf sets as input. cf. colour_tree_parts()."""
	colour_tree_parts(tree, constraints=constraints, **kw)
			
def colour_tree_with_leaf_groups(tree, leafgroups, **kw):
	"""paints arbitrary sets of leaves with a colour pallette; take a tree and a list of leaf sets as input. cf. colour_tree_parts()."""
	colour_tree_parts(tree, leafgroups=leafgroups, **kw)

## BioPython-using functions

def findSeqRecordIndexesFromSeqNames(aln, seqnames):
	if isinstance(seqnames, list): return [k for k,seq in enumerate(aln) if seq.id in seqnames]
	else:
		for k,seq in enumerate(aln):
			if seq.id in seqnames: return k

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

def parseChain(lnfchains, dold2newname={}, nfchainout=None, inchainfmt='nexus', outchainfmt='newick', verbose=False):
	"""parse a Nexus-format tree chain and re-write it in a Newick format; edit the tree on the fly, substituting tip labels or grafting trees on tips"""
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
