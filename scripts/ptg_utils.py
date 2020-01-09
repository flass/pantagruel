#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import tree2
from Bio import SeqIO, AlignIO, Align, Alphabet
from BCBio import GFF
from Bio.Phylo import BaseTree, NewickIO, NexusIO, _io as PhyloIO
from StringIO import StringIO
from random import randint
import gzip
import pipes, tempfile
import copy

supported_formats = {'newick': NewickIO, 'nexus': NexusIO}
dalnext = {'fasta':'fasta', 'nexus':'nex'}

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

def quantile(seq, P, ignoreNull=True):
	"""given seq a vector of values and P a vector of break probabilities, return the values of quantiles of seq at breaks
	
	assumes a continuity of the distribution as following a function linear by segment to approximate quantiles
	"""
	Q = []
	for p in P:
		l = [k for k in seq if ((k is not None) or (not ignoreNull))]
		l.sort()
		L = float(len(l))
		if not l: return None
		rankq = (L - 1)*p
		irankq = int(rankq)
		drankq = rankq - irankq
		if drankq==0:
			q = float(l[irankq])
		else:
			q = float(l[irankq])*(1-drankq) + float(l[irankq+1])*drankq
		Q.append(q)
	return tuple(Q)

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
	return [leaflab.rsplit(species_sep, 1)[0] for leaflab in leaflabs]

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
	force = kw.get('force', False)
	if ('constraints' in kw) and ('leafgroups' in kw):
		raise ValueError, "cannot pass both arguments 'constraints' and 'leafgroups' at the same time"
	leafsets = kw.get('constraints', kw.get('leafgroups'))
	cols = kw.get('palette', [colorWheel() for i in range(len(leafsets))])	
	if ('constraints' in kw):
		for k, cons in enumerate(reversed(leafsets)):
			# start from the last assigned node constraint, so that ealier nested clades will be painted over at a further iteration
			n = tree.mrca([dalias.get(ll, ll) for ll in cons], force=force)
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

def findCladeOutgroup(constraint, tree, didseq):
	anccons = tree.mrca(constraint, force=True) # constraint contain labels from the identical leaf set
	try:
		outgroup = anccons.go_brother()
		outgroupleaf = outgroup.get_leaf_labels()[0]
	except ValueError, e:
		sys.stderr.write( "Warning: in family %s, constraint %s leaf set makes a paraphyletic group.\nconstraint: %s\nanccons is root: %s\n"%(nffullali, cladename, constraint, str(anccons.is_root())) )
		for refidseq, redidseqs in didseq.iteritems():
			if (set(constraint) & set(redidseqs)):
				sw = "Inclusion of identical leaf set: {%s:%s}\n to constraint results in paraphyletic group.\n"%(refidseq, repr(redidseqs))
				sw+= "This indicate bad rooting of the gene tree, but can be ignored here.\nOutgroup for consraint %s will be set to 'None'.\n"%cladename
				sys.stderr.write( sw )
				outgroup = None
				outgroupleaf = 'None'
				break
		else:
			sys.stderr.write( "This indicate a definition of constraint clade that is inconsistent with the tree.\n" )
			raise ValueError, e
	return outgroupleaf

def openwithfilterpipe(filepath, mode, maskchars=None):
	"""opens file trough a character translating pipe.
	
	If maskchars evaluates to True, only 'r' or 'w' file opening modes are allowed due to limitation.
	in pipes.Templates.open(); specifying write-and-read mode 'w+' will result in write-only mode 'w'.
	"""
	if maskchars:
		filterpipe = pipes.Template()
		if len(maskchars[0])==1 and len(maskchars[1])==1:
			charfilter = "tr '%s' '%s'"%maskchars
		else:
			charfilter = "sed -e 's/%s/%s/g'"%maskchars
		filterpipe.append(charfilter, '--')
		f = filterpipe.open(filepath, mode.rstrip('+'))
	else:
		f = open(filepath, mode)
	return f

def treewrite(trees, filepath, treeformat, maskchars=None, branch_length_only=True, **kwargs):
	"""extended version of Phylo._io.write() allowing to translate of characters
	
	only does not support file handles (i.e. already opened files) as input, only file paths.
	"""
	if isinstance(trees, (BaseTree.Tree, BaseTree.Clade)):
		# Passed a single tree instead of an iterable -- that's OK
		trees = [trees]
	with openwithfilterpipe(filepath, 'w+', maskchars) as fp:
		n = getattr(supported_formats[treeformat], 'write')(trees, fp, branch_length_only=branch_length_only, **kwargs)
	return n

def treeappend(trees, filepath, treeformat, maskchars=None, branch_length_only=True, **kwargs):
	"""extended version of Phylo._io.write() allowing to use append file opening mode and to translate of characters
	
	only does not support file handles (i.e. already opened files) as input, only file paths.
	"""
	if isinstance(trees, (BaseTree.Tree, BaseTree.Clade)):
		# Passed a single tree instead of an iterable -- that's OK
		trees = [trees]
	if not maskchars:
		with openwithfilterpipe(filepath, 'a+') as fp:
			n = getattr(supported_formats[treeformat], 'write')(trees, fp, **kwargs)
	else:
		with open(filepath, 'a+') as fp:
			with tempfile.NamedTemporaryFile('w+') as t:
				with openwithfilterpipe(t.name, 'w', maskchars) as tfp:
					n = getattr(supported_formats[treeformat], 'write')(trees, tfp, branch_length_only=branch_length_only, **kwargs)
					tfp.flush()
				t.seek(0, 0)
				fp.write(t.read())
	return n

def treeparse(filepath, treeformat, maskchars=None, **kwargs):
	"""extended version of Phylo._io.parse() allowing to translate characters
	
	only does not support file handles (i.e. already opened files) as input, only file paths.
	"""
	with openwithfilterpipe(filepath, 'r', maskchars) as fp:
		for tree in getattr(supported_formats[treeformat], 'parse')(fp, **kwargs):
			yield tree

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

def parseChain(lnfchains, dold2newname={}, nfchainout=None, inchainfmt='nexus', outchainfmt='newick', maskchars=None, dirout='', verbose=False):
	"""parse a (Nexus-format) tree chain file and re-write it in a Newick format; edit the tree on the fly, substituting tip labels or grafting trees on tips
	
	A character filter is added as a file pipe when parsing the file, replacing any character of maskchars[0] with the counterpart in maskchars[1].
	The filter is reverted when writing the file. 
	With the filter maskchars=[('\([A-Z0-9]\)-', '\\1@'), ('\([A-Z0-9]\)@', '\\1-')] (using sed command), dash characters '-' PRECEDED BY A CAPITAL LETTER OR DIGIT 
	are translated into at symbols '@' on input, and '@'s are translated into '-'s on output; 
	this is to handle the fact that Bio.Nexus tree parser does not support dashes in the taxon labels (see https://github.com/biopython/biopython/issues/1022).
	!!! CAUTION: any '@' in the original label will thus be turned into a '-' in the final file output. Please avoid '@'s in tree taxon labels. !!!
	!!! CAUTION 2: maskchars=('-', '@') leads to wrong behaviour if hyphen float numbers noted Xe-XXX are substituted as leads to branch lengths to be read as names.
	"""
	if maskchars:
		if isinstance(maskchars, list):
			invmaskchars = maskchars[1]
			maskchars = maskchars[0]
		else:
			invmaskchars = (maskchars[1], maskchars[0])
	else:
		invmaskchars = None
	if verbose:
		print "parseChain('%s')"%repr(lnfchains)
		if maskchars: print "'%s' character set will be substituted with '%s' on input and back on output"%maskchars
	ntree = 0
	buffsize = 50
	treebuffer = []
	dhandles = {}
	# prepare handles for multiple input files
	for k, nfchain in enumerate(lnfchains):
		# by default add a pipe filter to deal with lack of support for dash characters in taxon labels in BioPython Nexus tree parser
		dhandles[k] = treeparse(nfchain, inchainfmt, maskchars=maskchars)
	# pepare output file
	if nfchainout: nfout = nfchainout
	else: nfout = n2outChainFileName(lnfchains[0], dirout, inchainfmt, outchainfmt)
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
								raise ValueError, "replacement value for a leaf must be either a string (to edit leaf label) or a tree object instance of classes tree2.Node or Bio.Phylo.TreeElement (or derivates)"
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
					treewrite(treebuffer, nfout, outchainfmt, maskchars=invmaskchars)
				else:
					treeappend(treebuffer, nfout, outchainfmt, maskchars=invmaskchars)
				treebuffer = []
	if treebuffer: treeappend(treebuffer, nfout, outchainfmt, maskchars=invmaskchars)
	print '%s%s ...done (%d trees)'%(('\n' if verbose else ''), nfout, ntree)
	return None

def replaceInSingleTree(nfgt, dold2newname={}, nfgtout=None, ingtfmt='newick', outgtfmt='newick', maskchars=None, verbose=False, dirout='', mode='tree2.Node'):
	"""(single-tree version of parseChain) parse a (Newick-format) tree file and re-write it in a Newick format; edit the tree on the fly, substituting tip labels or grafting trees on tips"""
	
	### alternative codes
	def replaceInSingleTree_withBioPhylo(nfgt, dold2newname, nfout, ingtfmt, outgtfmt, maskchars, verbose):
		"""Using Bio.Phylo - direct equivalent of parseChain
		
		complex code with no real point when not parsing large ammount of trees from each input file
		
		A character filter is added as a file pipe when parsing the file, replacing any character of maskchars[0] with the counterpart in maskchars[1].
		The filter is reverted when writing the file. 
		With the filter maskchars=[('\([A-Z0-9]\)-', '\\1@'), ('\([A-Z0-9]\)@', '\\1-')] (using sed command), dash characters '-' PRECEDED BY A CAPITAL LETTER OR DIGIT 
		are translated into at symbols '@' on input, and '@'s are translated into '-'s on output; 
		this is to handle the fact that Bio.Nexus tree parser does not support dashes in the taxon labels (see https://github.com/biopython/biopython/issues/1022).
		!!! CAUTION: any '@' in the original label will thus be turned into a '-' in the final file output. Please avoid '@'s in tree taxon labels. !!!
		!!! CAUTION 2: maskchars=('-', '@') leads to wrong behaviour if hyphen float numbers noted Xe-XXX are substituted as leads to branch lengths to be read as names.
		"""
		if maskchars:
			if isinstance(maskchars, list):
				invmaskchars = maskchars[1]
				maskchars = maskchars[0]
			else:
				invmaskchars = (maskchars[1], maskchars[0])
		else:
			invmaskchars = None
		if verbose:
			print "replaceInSingleTree('%s')"%repr(lnfgt)
			if maskchars: print "'%s' character set will be substituted with '%s' on input and back on output"%maskchars
		dhandles = {}
		# parse input file
		treehandle = treeparse(nfgt, ingtfmt, maskchars=maskchars)
		tree = treehandle.next()
		if dold2newname:
			for tip in tree.get_terminals():
				nutipname = dold2newname.get(tip.name)
				if nutipname:
					if isinstance(nutipname, str):
						if verbose: print "replace leaf label '%s' with label %s"%(tip.name, repr(nutipname))
						tip.name = nutipname
					else:
						if isinstance(nutipname, tree2.Node):
							nusubtree = tree2toBioPhylo(nutipname)
						elif isinstance(nutipname, BaseTree.TreeElement):
							nusubtree = nutipname
						else:
							raise ValueError, "replacement value for a leaf must be either a string (to edit leaf label) or a tree object instance of classes tree2.Node or Bio.Phylo.TreeElement (or derivates)"
						if verbose: print "replace leaf labelled '%s' with clade %s"%(tip.name, repr(nutipname))
						subtreelen = nusubtree.total_branch_length()/nusubtree.count_terminals()
						if subtreelen:
							# substract the subtree length to its branch length
							tip.branch_length = max(0.0, tip.branch_length - subtreelen)
						# attach subtree to tree
						tip.clades = nusubtree.root.clades
						tip.name = ''
		treewrite(tree, nfout, outgtfmt, maskchars=invmaskchars)
		if verbose: print '\n%s ...done'%(nfout)
		return None
	
	def replaceInSingleTree_withtree2Node(nfgt, dold2newname, nfgtout, ingtfmt, outgtfmt, verbose):
		"""Using tree2.Node - no support for replacement clades in Bio.Phylo.TreeElement format
		
		less object format conversions and more straightforward code so more efficient when only a single trees is read from each input file
		"""
		tree = getattr(tree2, 'read_%s'%ingtfmt)(nfgt)
		if dold2newname:
			for tip in tree.get_leaves():
				nutipname = dold2newname.get(tip.label())
				if nutipname:
					if isinstance(nutipname, str):
						if verbose: print "replace leaf label '%s' with label %s"%(tip.label(), repr(nutipname))
						tip.edit_label(nutipname)
					else:						
						if verbose: print "replace leaf labelled '%s' with clade %s:"%(tip.label(), repr(nutipname))
						if not isinstance(nutipname, tree2.Node):
							raise ValueError, "replacement value for a leaf must be either a string (to edit leaf label) or a tree object instance of classes tree2.Node or Bio.Phylo.TreeElement (or derivates)"
						nusubtree = nutipname.deepcopybelow()
						nusubtree.edit_label('')
						tipfat = tip.go_father()
						if verbose:
							print tipfat
							print '  ->'
						subtreelen = nusubtree.treelength()/nusubtree.nb_leaves()
						# substract the subtree length to its branch length
						newlen = max(0.0, tip.lg() - subtreelen)
						newboot = tip.bs()
						# remove tip node from tree
						tipfat.unlink_child(tip)
						# attach subtree node to tree
						tipfat.link_child(nusubtree, newlen=newlen, newboot=newboot)
						if verbose: print tipfat
		getattr(tree2, 'write_%s'%outgtfmt)(tree, nfout)
		if verbose: print '\n%s ...done'%(nfout)
		return None

	### run
	# pepare output file
	if nfgtout: nfout = nfgtout
	else: nfout = n2outChainFileName(nfgt, dirout, ingtfmt, outgtfmt)
	if mode=='Bio.Phylo':
		return replaceInSingleTree_withBioPhylo(nfgt, dold2newname, nfout, ingtfmt, outgtfmt, maskchars, verbose)
	elif mode=='tree2.Node':
		return replaceInSingleTree_withtree2Node(nfgt, dold2newname, nfout, ingtfmt, outgtfmt, verbose)

def labsFromReplacementLabOrSubtree(newlaborst):
	if isinstance(newlaborst, str):
		newlabs = [newlaborst]
	elif isinstance(newlaborst, tree2.Node):
		newlabs = newlaborst.get_leaf_labels()
	elif isinstance(newlaborst, BaseTree.TreeElement):
		newlabs = [tip.name for tip in newlaborst.get_terminals()]
	else:
		raise ValueError: "unsupported format for replacement label or subtree: %s"%repr(newlaborst)
	return newlabs

def duplicateSeqsInAln(nfcolaln, dold2newname, nfoutreplaln=None, inalnfmt='nexus', outalnfmt='fasta'):
	colaln = AlignIO.read(nfcolaln, format=inalnfmt, alphabet=Alphabet.generic_dna)
	replaln = Align.MultipleSeqAlignment([], alphabet=Alphabet.generic_dna)
	rmseqrowid = []
			
	for cladename, newlaborst in dold2newname.iteritems():
		# collect collpased clade representative sequences
		ccreprseqalnrowid = findSeqRecordIndexesFromSeqNames(aln, cladename)
		rmseqrowid.append(ccreprseqalnrowid)
		ccreprseq = colaln[ccreprseqalnrowid]
		newlabs = labsFromReplacementLabOrSubtree(newlaborst)
		for newlab in newlabs:
			newseq = copy.copy(ccreprseq) # shallow copy is enough as 'id' attribute is a str i.e. non referenceable object
			newseq.id = newlab
			replaln.append(newseq)
	
	for rowseqid in range(len(colaln):
		if not rowseqid in rmseqrowid:
			replaln.append(colaln[rowseqid])
	
	if nfoutreplaln: AlignIO.write(replaln, nfoutreplaln, outalnfmt)
	return replaln

def seqrecordsFromGBFF(nfgbff):
	if nfgbff.endswith('.gz'):
		fgbff = gzip.open(nfgbff, 'rb')
	else:
		fgbff = open(nfgbff, 'r')
	genome = SeqIO.read(fgbff, 'genbank')
	fgbff.close()
	return genome

def seqrecordsFromGFFandGenomicFasta(nfgff, nffastain):
	if nfgff.endswith('.gz'):
		fgff = gzip.open(nfgff, 'rb')
	else:
		fgff = open(nfgff, 'r')
	if nffastain.endswith('.gz'):
		ffastain = gzip.open(nffastain, 'rb')
	else:
		ffastain = open(nffastain, 'r')
	seqdict = SeqIO.to_dict(SeqIO.parse(ffastain, "fasta", alphabet=Alphabet.generic_dna))
	genome = list(GFF.parse(fgff, seqdict))
	fgff.close()
	ffastain.close()
	return genome

def extractCDSFasta(seqrecord, feature, ffastaout, ncds):
	ncds += 1
	recid = seqrecord.id
	cdsseq = feature.location.extract(seqrecord).seq
	protid = feature.qualifiers.get('protein_id', [None])[0]
	qualifs = ' '.join("[%s=%s]"%(str(k), str(v[0])) for k,v in feature.qualifiers.iteritems())
	if protid:
		ffastaout.write(">lcl|%s_cds_%s_%d %s\n%s\n" % (recid, protid, ncds, qualifs, cdsseq) )
	else:
		ffastaout.write(">lcl|%s_cds_%d %s %s\n%s\n" % (recid, ncds, qualifs, cdsseq) )
	return ncds

def extractCDSFastaFromSeqrecords(genome, nffastaout):
	if nffastaout.endswith('.gz'):
		ffastaout = gzip.open(nffastaout, 'wb')
	else:
		ffastaout = open(nffastaout, 'w')
	ncds = 0
	for seqrecord in genome:
		for feature in seqrecord.features:
			if feature.type == "CDS":
				ncds = extractCDSFasta(seqrecord, feature, ffastaout, ncds)
			if hasattr(feature, 'sub_features'):
				# specific to records from GFF.parse
				loctagqual = feature.qualifiers.get('locus_tag', [''])
				for subfeature in feature.sub_features:
					if subfeature.type == "CDS":
						subfeature.qualifiers['locus_tag'] = loctagqual
						ncds = extractCDSFasta(seqrecord, subfeature, ffastaout, ncds)
	ffastaout.close()

def extractCDSFastaFromGFFandGenomicFasta(nfgff, nffastain, nffastaout):
	genome = seqrecordsFromGFFandGenomicFasta(nfgff, nffastain)
	extractCDSFastaFromSeqrecords(genome, nffastaout)
	
def extractCDSFastaFromGBFF(nfgbff, nffastaout):
	genome = seqrecordsFromGBFF(nfgbff)
	extractCDSFastaFromSeqrecords(genome, nffastaout)

#### Database connection functions

# allow seemless transition between db engines
def connectpostgresdb(dbname, **kw):
	psycopg2 = __import__('psycopg2')
	return psycopg2.connect(dbname=dbname, **kw)

def connectsqlitedb(dbname):
	sqlite3 = __import__('sqlite3')
	return sqlite3.connect(dbname)

def get_dbconnection(dbname, dbengine):
	if (dbengine.lower() in ['postgres', 'postgresql', 'psql', 'pg']):
		dbcon = connectpostgresdb(dbname)
		dbtype = 'postgres'
		valtoken='%s'
		dbcur = dbcon.cursor()
		dbcur.execute("set search_path = genome, phylogeny;")
	elif (dbengine.lower() in ['sqlite', 'sqlite3']):
		dbcon = connectsqlitedb(dbname)
		dbcur = dbcon.cursor()
		dbtype = 'sqlite'
		valtoken='?'
	else:
		raise ValueError,  "wrong DB type provided: '%s'; select one of dbengine={'postgres[ql]'|'sqlite[3]'}"%dbengine
	return (dbcon, dbcur, dbtype, valtoken)

