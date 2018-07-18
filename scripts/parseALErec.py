#!/usr/bin/python
import tree2
import sys, os, getopt
import re
import copy

eventTypes = 'DTLS'

def getOriSpeciesFromEventLab(eventlab, sgsep='_'):
	# split at DT location separator '@', then possibly at T don/rec seprator '->', and finally shorten the species label if node is a leaf
	elab = eventlab.split('@')[1] if '@' in eventlab else eventlab
	return elab.split('->')[0].split(sgsep)[0]

def parseALERecFile(nfrec, reftreelen=None, restrictclade=None, skipEventFreq=False):
	line = ''
	lrecgt = []
	restrictlabs = []
	frec = open(nfrec, 'r')
	while not line.startswith('S:\t'):
		line = frec.readline()
	# extract node labels from reconciled species tree
	spetree = tree2.AnnotatedNode(nwk=line.strip('\n').split('\t')[1], namesAsNum=True)
	spetree.complete_node_ids()
	if reftreelen:
		if not spetree.hasSameTopology(reftreelen): raise IndexError, "reference tree from $2 has not the same topology as that extracted from reconciliation output file $1"
		for node in spetree:
			# extract branch length from topologically identical tree from $2
			matchclade = reftreelen.map_to_node(node.get_leaf_labels())
			node.set_lg(matchclade.lg())
	if restrictclade:
		for restrictnodelab in restrictclade.split(','):
			restrictlabs += spetree[restrictnodelab].get_children_labels()
		subspetree = spetree.restrictToLeaves(restrictlabs, force=True)
	else:
		subspetree = spetree
	while not line.endswith('reconciled G-s:\n'):
		line = frec.readline()
	for i in range(2): line = frec.readline() # skips 2 lines
	# extract reconciled gene tree(s)
	recgtlines = []
	while not line.startswith('#'):
		recgtlines.append(line)
		rectree = tree2.AnnotatedNode(nwk=line.strip('\n'), namesAsNum=True)
		rectree.complete_node_ids()
		lrecgt.append(rectree)
		line = frec.readline()
	dnodeevt = {}
	if not skipEventFreq:
		for i in range(3): line = frec.readline() # skips 3 lines
		# extract node-wise event frequency / copy number info
		for line in frec:
			if line=='\n': continue
			lsp = line.strip('\n').split('\t')
			dnodeevt[lsp[1]] = [float(s) for s in lsp[2:]]
	frec.close()
	return [spetree, subspetree, lrecgt, recgtlines, restrictlabs, dnodeevt]

def parseUndatedRecGeneTree(recgt, spet, dexactevt={}, recgtsample=None, nsample=1, sgsep='_', restrictlabs=[], fillDTLSdict=True, recordEvTypes='DTL', excludeTaggedLeaves=None, excludeTaggedSubtrees=None):
	"""extract list of events from one reconciled gene tree as found in output of ALEml_undated (Szolosi et al., 2013; https://github.com/ssolo/ALE)
	
	This function returns two objects:
	
	- a dict object 'dlevt', containing the (possibly redundant) list of unique events having occurred in the tree,
	sorted by event type, which can be 'D', 'T', 'L' or 'S', for duplications, horizontal transfers, losses or speciations.
	This record is a crude aggregate of what happened in the whole gene family, present for back-compatiility purposes.
	The behaviour of filling up this record can be turned off to gain time and memory with fillDTLSdict=False.
	
	- a dict object 'dnodeallevt', containing the list all events in one reconciliation scenario
	(a single scenario anomg the sample), sorted by gene tree node id. 
	Event are represented as tuples of the following form: (X, [don, ] rec, freq),
	where X is the event type
	'rec' and optionally 'don' (for Ts only) are species tree node labels where the event where inferred,
	and 'freq' the event frequency of this event in the whole sample.
	
	Frequency of events is searched in the WHOLE reconciled gene tree sample 
	provided as the list of pseudo-newick strings, using exact string matching (only for DT).
	In case of the same pattern of event occurring repeatedly in the same tree, e.g. same T@DON->REC in two paralogous lineages 
	(can likely happen with tandem duplicates...), the count will reflect the sum of all such events. 
	Records of event frequencies in both 'dlevt' and 'dnodeallevt' are thus NOT differentiated by lineage.
	"""
	def parseRecGTNode(node, dlevt, dnodeallevt):
		nodelab = node.label()
		if not nodelab: raise ValueError, "unannotated node:\n%s"%str(node)
		if node.is_leaf() and (excludeTaggedLeaves in nodelab):
			# the species assignment of this leaf is not certain (inferred)
			# and thus events leading directly to this leaf are not to be trusted and should be skipped
			return
		nodeid = node.nodeid()
		# line of events to be read left-to-right backward in time
		lineage = nodelab.split('.')
		for i in range(1, len(lineage)):
			eventlab = lineage[i]
			# identify next (i.e. forward in time) event on the lineage
			preveventlab = lineage[i-1] 
			if eventlab.startswith('D@'):
				if ('D' in recordEvTypes):
					# duplication event
					dup = eventlab.split('D@')[1]
					evtup = ('D', dup)
					fevt = dexactevt.setdefault(evtup, float(recgtsample.count(eventlab))/nsample)
					if restrictlabs and not (dup in restrictlabs): continue
					if fillDTLSdict: dlevt['D'].append(dup)
					dnodeallevt.setdefault(nodeid, []).append(evtup)
			elif eventlab.startswith('T@'):
				if ('T' in recordEvTypes):
					# transfer event
					translab = eventlab.split('T@')[1]
					don, rec = translab.split('->')
					evtup = ('T', don, rec)
					fevt = dexactevt.setdefault(evtup, float(recgtsample.count(eventlab))/nsample)
					if restrictlabs and not ((don in restrictlabs) and (rec in restrictlabs)): continue
					if fillDTLSdict: dlevt['T'].append((don, rec))
					dnodeallevt.setdefault(nodeid, []).append(evtup)
			else:
				# just a species tree node label
				if ('S' in recordEvTypes):
					spe = getOriSpeciesFromEventLab(eventlab, sgsep=sgsep)
					evtup = ('S', spe)
					#~ evtpat = "(\)\.|[\(,])%s"%eventlab	# pattern captures the event at an internal node or a leaf 
					evtpat = "([\.\(,]%s)"%eventlab	# pattern captures the event at an internal node or a leaf 
					#~ print "count occurence in sample of pattern:", evtpat
					fevt = dexactevt.setdefault(evtup, float(len(re.search(evtpat, recgtsample).groups()))/nsample)
					if fillDTLSdict: dlevt['S'].append(spe)
					dnodeallevt.setdefault(nodeid, []).append(evtup)
				if preveventlab!='':
					if ('L' in recordEvTypes):
						# speciation-loss event
						# speciation occurs of the named node but loss acutally occurs in its descendant (the other than that below/preceding on the lineage)
						ploss = spet[eventlab]
						closslabs = ploss.children_labels()
						closslabs.remove(getOriSpeciesFromEventLab(preveventlab, sgsep=sgsep))
						if len(closslabs)>1: raise IndexError, "non binary species tree at node %s (children: %s)"%(lineage[-1], repr(ploss.get_children_labels()))
						los = closslabs[0]
						evtup = ('L', los)
						fevt = dexactevt.setdefault(evtup, float(re.search("([^\)]\.%s)"%eventlab, recgtsample).groups)/nsample)
						if restrictlabs and not (los in restrictlabs): continue
						if fillDTLSdict: dlevt['L'].append(los)
						dnodeallevt.setdefault(nodeid, []).append(evtup)
				#~ else:
					#~ # a simple speciation event ; already delt with
					#~ pass
		if excludeTaggedSubtrees:
			alllabext = [lab.split('_', 1)[1] for lab in node.iter_leaf_labels()]
			if (len(set(alllabext)) == 1) and (excludeTaggedSubtrees in alllabext[0]):
				# the subtree below this node is an artificial addition to the reconciled gene tree and should be skipped
				return
		for child in node.get_children():
			# use recursion to be able to exclude subtrees
			parseRecGTNode(child, dlevt, dnodeallevt)
		return
	
	dnodeallevt = {}
	dlevt = {e:[] for e in eventTypes}
	parseRecGTNode(recgt, dlevt, dnodeallevt) # recursive function call
	return dlevt, dnodeallevt

def parseRecGeneTree(recgt, spet, ALEmodel='undated', **kw):
	if ALEmodel=='undated':
		return parseUndatedRecGeneTree(recgt, spet, **kw)
	elif ALEmodel=='dated':
		return parseDatedRecGeneTree(recgt, spet, **kw)

#### dated model gene tree syntax parsing
def findNextDatedEventSlice(nodelab): #, receptionOK=False
	nextev1 = nodelab.find('@')
	if nextev1 > 0:
		if nodelab[nextev1-1] in {'T', 'D', 'b'}:
			nextev1 -= 1		# account for the preceding D or T
			if nodelab[nextev1-1]=='b':
				if nodelab[nextev1-2]=='T':
					nextev1 -= 2	# account for the preceding Tb
				else:
					raise IndexError
		#~ elif not receptionOK:
			#~ raise IndexError
	nextev2 = nodelab.find('.')
	if (nextev1 < 0) and (nextev2 < 0):
		return None
	elif (nextev1 >= 0) and ((nextev2 < 0) or (nextev1 < nextev2)):
		return nextev1
	elif (nextev2 >= 0) and ((nextev1 < 0) or (nextev2 < nextev1)):
		return nextev2
	else:
		raise IndexError

def splitSingleUndatedEvent(nodelab, isleaf=False):
	"""yet to implement; just need to check really as in principle splitting at periods is enough"""
	raise NotImplementedError
	lineage = nodelab.split('.') # bring in code from parseUndatedRecGeneTree()
	return lineage

def extractLabelfromDatedEventLeaf(nodelab):
	"""extract leaf label prefix from string (possible followed by an event chain)"""
	nextev = findNextDatedEventSlice(nodelab) #, receptionOK=True)
	leaflab = nodelab[:nextev] # the whole string if nextev is None
	return leaflab

def splitSingleDatedEvent(nodelab, isleaf=False, verbose=False, **kw):
	"""from a gene tree branch event chain, return iterator that yields the type and location of every event
	
	if isleaf=True, assumes the branch led to a leaf, and the leaf label is yielded as the last item.
	"""
	def process_event(s):
		# process possible event separators in order
		if s.startswith('.T@'):
			# immeadiate transfer-loss (should be only present on the beginning of the event chain)
			# emission of transfer by donor followed by loss of resident copy; just changes species identity of the gene tree branch (from donor to recipent species)
			# will be followed by transfer reception event
			thisev = 3
			evtype = 'TdL'
		elif s.startswith('.'):
			# simple speciation
			evtype = 'S'
			thisev = 1
		elif s.startswith('@'):
			# reception of transfer by recipient = gain
			thisev = 1
			evtype = 'Tr'
		elif s.startswith('Tb@'):
			# reception of transfer by recipient (case where diversification occurs out of ref tre and several copies come back in) = gain
			thisev = 3
			evtype = 'Tr'
		elif s.startswith('T@'):
			# emission of transfer by donor; not relevant to gene content
			thisev = 2
			evtype = 'Td'
		elif s.startswith('D@'):
			# duplication = gain
			thisev = 2
			evtype = 'D'
		else:
			raise ValueError
		nextev = findNextDatedEventSlice(s[thisev:]) #, receptionOK=reception)
		if nextev is None:
			# no other event
			evlocdate = s[thisev:]
			s = None
		else:
			evlocdate = s[thisev:nextev+thisev]
			s = s[nextev+thisev:]
			# if event is not last it is followed by a loss: loss in one of the daughter lineage for speciation; loss of resident copy for transfer
			# just changes species identity of the gene tree branch (from parent to daughter species)
			if evtype=='S': evtype += 'L'
		if evtype in ['S', 'SL']:
			evloc = evlocdate
			evdate = evlocdate
		else:
			evdate, evloc = evlocdate.split('|')
		return s, evtype, evloc, int(evdate)
	
	# initiate
	if isleaf:
		leaflab = extractLabelfromDatedEventLeaf(nodelab)
		# remove leaf label prefix before parsing (possible) trailing event chain
		s = nodelab[len(leaflab):]
		if verbose>1: print 'leaflab: %s;'%leaflab,
	else:
		s = nodelab
	# parse string describing the event chain
	if verbose>1: print 'events:',
	while s:
		s, evtype, evloc, evdate = process_event(s)
		yield (evtype, evloc, evdate)
		if verbose>1: print (evtype, evloc, evdate),
	if isleaf:
		# add the final speciation
		evtype = 'S'
		evdate = 0
		evloc = leaflab
		yield (evtype, evloc, evdate)
		if verbose>1: print 'leaflab: %s;'%leaflab,
		# last yield
		yield leaflab
	if verbose>1: print ''

def splitEventChain(nodelab, isleaf=False, ALEmodel='dated', **kw):
	"""from the event chain string, return the list of event tuples (type, date, location) forward in time, but for last element that is the leaf label (or None if was not a leaf)"""
	verbose = kw.get('verbose')
	if ALEmodel=='undated':
		lineage = [event for event in splitSingleUndatedEvent(nodelab, isleaf=isleaf, **kw)]
	elif ALEmodel=='dated':
		lineage = [event for event in splitSingleDatedEvent(nodelab, isleaf=isleaf, **kw)]
		# line of events is read left-to-right FORWARD in time; reverse it for a BACKWARD result
		lineage.reverse()
	if isleaf:
		# remove the first-yielded element from event list that's the leaf label
		leaflab = lineage[0]
		lineage = lineage[1:]
	else:
		leaflab = None
	if verbose: print 'lineage:', lineage, ("leaflab: %s"%leaflab if leaflab else "")
	return (lineage, leaflab)

def _prune_orthologs_bottom_up(node, ALEmodel='dated', **kw):
	"""'last-gain' (strict) definition of ortholous groups: only those genes related by a line of speciation events (no transfer, no duplication) are orthologs"""
	# initiate at root
	unclassified = set([])
	orthologGroups = kw.get('orthologGroups', [])
	dlabs = kw.get('dlabs', {})
	verbose = kw.get('verbose')
	nodelab = node.label()
	if not nodelab: raise ValueError, "unannotated node:\n%s"%str(node)
	lineage, leaflab = splitEventChain(nodelab, isleaf=node.is_leaf(), **kw)
	if node.is_leaf():
		unclassified.add(leaflab) # where unclasified set is filled up
		dlabs[nodelab] = leaflab
	else:
		# post-order traversal recursion:
		for child in node.children:
			# first explore children
			#~ orthologGroups, uncl, dlabs = _prune_orthologs_bottom_up(child, orthologGroups=orthologGroups, dlabs=dlabs, **kw)
			for arg in ('orthologGroups', 'dlabs'): kw[arg] = eval(arg)
			orthologGroups, uncl, dlabs = _prune_orthologs_bottom_up(child, **kw)
			unclassified |= uncl
	# then explore node itself
	if not unclassified:
		# all leaves below have been already classified in orthologous group; stop here.
		return orthologGroups, unclassified, dlabs
	# list of events goes BACKWARD in time when read left-to-right
	for event in lineage:
		evtype, evloc, evdate = event
		if evtype in {'Tr', 'T'}:
			# gene was last gained here by a species ancestor:
			# what was not yet classifed in this subtree is an orthologous group
			ortho = tuple(sorted(unclassified))
			orthologGroups.append(ortho)
			unclassified = set([])
			if verbose: print 'T-OG!'
			break # for event loop
		elif evtype=='D':
			# gene was duplicated here by species ancestor:
			# what was not yet classifed in the child subtrees
			# of this subtree are each an orthologous group
			for child in node.children:
				subuncl = set([dlabs[leaflabevchain] for leaflabevchain in child.iter_leaf_labels()]) & unclassified
				ortho = tuple(sorted(subuncl))
				orthologGroups.append(ortho)
			unclassified = set([])
			if verbose: print 'D-OG! D-OG!'
			break # for event loop
	#~ if verbose: print "unclassified:", unclassified
	if node.is_root():
		# any remaining unclassified leaf is to be alloacted to a backbone orthologous group
		# (= paraphyletic gene lineage unaffected by transfers since the origin of the family)
			ortho = tuple(sorted(unclassified))
			orthologGroups.append(ortho)
			if verbose: print 'backbone OG!'
	return orthologGroups, unclassified, dlabs

def getSpe(lab, sp0, sp1):
	return lab.split(sp0, sp1)[0]

def getExtraNumerarySpe(lspe, ancclade=[]):
	sspe = set(lspe)
	if ancclade: baselist = ancclade
	else: baselist = list(sspe)
	return reduce(lambda x,y: x+y, [[spe]*(lspe.count(spe) - baselist.count(spe)) for spe in sspe], [])

def _prune_nested_candidate_orthologs(leaflabs, lspe, extraspe, candidateOGs, refspetree=None, anclade=None, **kw):
	""""""
	# first filter out the last-gain defined OGs that do not map completely under this clade 
	# this is to remove a potential "backbone" OG = extremely paraphyletic group 
	# that is the remainder of the tree after the removal of all others last-gain OGs
	sp0, sp1 = kw.get('splitparam', ('_',1))
	verbose = kw.get('verbose')
	orthologGroups = []
	sleaflabs = set(leaflabs)
	cOGs = [set(cOG) for cOG in candidateOGs if set(cOG) <= sleaflabs]
	# score the last-gain defined OGs based on their capacity to remove the extra species or copies from the leaf set
	sextraleaflabs = set([leaflab for leaflab in leaflabs if (getSpe(leaflab, sp0, sp1) in extraspe)])
	cOGs.sort(key=lambda x: len(sextraleaflabs & x), reverse=True) # sort first the cOG with the hihest score of leaf set overlap
	clspe = copy.copy(lspe)
	for k, cOG in enumerate(cOGs):
		if verbose: print 'try substracting last-gain OG: %s'%repr(tuple(cOG))
		for ll in cOG:
			clspe.remove(getSpe(ll, sp0, sp1))
		if refspetree: extraspe = getExtraNumerarySpe(clspe, ancclade)
		else: extraspe = getExtraNumerarySpe(lspe)
		if not extraspe:
			# first add the retained last-gain OGs to the final list of OGs
			orthologGroups += [tuple(sorted(cOG)) for cOG in cOGs[:(k+1)]]
			# then substract them from the current clade
			for cOG in cOGs[:(k+1)]:
				sleaflabs -= cOG
			ortho = tuple(sorted(sleaflabs))
			orthologGroups.append(ortho)
			if verbose: print 'U-OG! - '+' '.join(['G-OG!']*(k+1))
			rcOGs = [cOG for cOG in candidateOGs if (not cOG in orthologGroups)]
			break # for cOG loop
	else:
		rcOGs = candidateOGs
	return orthologGroups, rcOGs

def _prune_orthologs_top_down(node, **kw):
	"""from a ALE-reconciled gene tree, return a list of orthologous groups and the dict of leaf labels to the actual gene sequence name (removing the trailing annotation of event chain)
	
	Orthologous groups (OGs) are defined as the gene tree clades containing at most one representative of each species.  
	This recognizes unicopy clades as OGs, even those where transfer events occurred, as long as they did not induce 
	change of copy number in represented pecies (when considering the sum of all events in the subtree), i.e. that 
	transfers were gene conversion events.
	
	if 'refspetree' is specified:
	    another constraint is added that those represented species must all belong to the species clade below the ancestor 
	    to which a gene gain is mapped on the gene tree branch (last recipient of the gene if several gain events occurred).
	if 'candidateOGs' is specified:
	    this argument allow to pass on other putative OGs (typically of smaller span), which the function attempts to prune 
	    from clades that are not unicopy so to obtain an unicopy group. This is useful to implement mixed criterion, where
	    for instance late gains (by transfer from a distant donor or duplication) can be removed to capture the orthologous 
	    backbone of a clade.
	"""
	# # # # initiate paramaters
	orthologGroups = kw.get('orthologGroups', [])
	dlabs = kw.get('dlabs', {})
	sp0, sp1 = kw.get('splitparam', ('_',1))
	refspetree = kw.get('refspetree')
	candidateOGs = kw.get('candidateOGs')
	verbose = kw.get('verbose')
	if verbose and candidateOGs: print 'using mixed criterion'
	if not dlabs:	
		# first establish the dictionary of actual leaf lables, without the trailing event chain
		for leaflabevchain in node.iter_leaf_labels():
			dlabs[leaflabevchain] = extractLabelfromDatedEventLeaf(leaflabevchain)
	
	# # # # test at current node
	nodelab = node.label()
	if not nodelab: raise ValueError, "unannotated node:\n%s"%str(node)
	print repr(node)+';', ('leaf' if node.is_leaf() else 'internal')+';',
	leaflabs = [dlabs[leaflabevchain] for leaflabevchain in node.iter_leaf_labels()]
	lspe = [getSpe(leaflab, sp0, sp1) for leaflab in leaflabs]
	extraspe = True # initiate
	if not refspetree:
		if verbose: print 'evaluate orthology under gene tree node', nodelab
		extraspe = getExtraNumerarySpe(lspe)
		if not extraspe:
			ortho = tuple(sorted(leaflabs))
			orthologGroups.append(ortho)
			if verbose: print 'U-OG!'
		elif candidateOGs:
			if verbose: print 'extraspe:', extraspe
			OGs, rcOGs = _prune_nested_candidate_orthologs(leaflabs, lspe, extraspe, candidateOGs, **kw)
			if OGs:
				orthologGroups += OGs
				candidateOGs = rcOGs
				extraspe = []
	else:
		lineage, tiplab = splitEventChain(nodelab, isleaf=node.is_leaf(), **kw)
		# list of events goes BACKWARD in time when read left-to-right
		for event in lineage:
			# latest event primes
			evtype, evloc, evdate = event
			if evtype in {'Td', 'TdL'}: 
				# emission of transfer: irrelevant to gene content of the clade below this ancestor
				continue # the for event loop
			speanc = refspetree[evloc]
			assert isinstance(speanc, tree2.Node)
			if verbose: print 'evaluate orthology under gene tree event', evtype, 'at ancestor', evloc
			ancclade = speanc.get_leaf_labels()
			extraspe = getExtraNumerarySpe(lspe, ancclade)
			if not extraspe:
				ortho = tuple(sorted(leaflabs))
				orthologGroups.append(ortho)
				if verbose: print 'U-OG!'
				break # for event loop
			elif candidateOGs:
				if verbose: print 'extraspe:', extraspe
				OGs, rcOGs = _prune_nested_candidate_orthologs(leaflabs, lspe, extraspe, candidateOGs, refspetree=refspetree, anclade=anclade, **kw)
				if OGs:
					orthologGroups += OGs
					candidateOGs = rcOGs
					extraspe = []
					break
	if extraspe:
		# if did not find the clade to be an orthologous group, recurse down the tree:
		for child in node.children:
			for arg in ('orthologGroups', 'dlabs', 'candidateOGs'): kw[arg] = eval(arg)
			orthologGroups, dlabs = _prune_orthologs_top_down(child, **kw)
	return orthologGroups, dlabs

def getOrthologues(recgt, ALEmodel='dated', method='last-gain', **kw):
	"""prune reconciled gene tree at transfer or duplication branches to identify orthologous groups (OGs) of genes. 
	
	mandatory argument: 'method'
	several methods are possible:
	'last-gain' (alias: 'strict'):
	    bottom-up approach (tip-to-root) where only those leaves related by a 
	    line of speciation events (no transfer, no duplication) are orthologs; 
	    the search is achieved by 
	'unicopy':
	    top-down approach (root-to-tips) where orthologous groups are defined as 
	    the gene tree clades containing only unique representatives of species.
	    This allows to recognize as orthologous groups clades where 
	    transfer events within the clade that do not induce change of 
	    copy number (when considering the sum of all events in the subtree), 
	    i.e. gene conversion events
	'mixed': (heuristic ?) approach mixing the two previous: 
	    when exploring the tree with the top-own inclusve approach, presence of
	    excess leaves causing rejection of the OG clade status can be salvaged by 
	    pruning of last-gain OG detected in the bottom-up approach may reveal 
	    bona-fide unicopy OGs (even though not clades because of excluded 
	    secondary gains, treated as other OGs).
	
	other optional keyword arguments:
	'refspetree':
	    tree object of class tree2.Node (or a derived class). This tree must be annotated with 
	    the same internal node label as used in the reconciliation (as from the header of ALE 
	    reconciliation file)
	    Useful only when using the 'unicopy' or 'mixed' methods, the reference species tree 
	    allow to constrain the scope of unicopy clades: represented species 
	    must all belong to the species clade below the ancestor to which a gene gain is mapped 
	    on the gene tree branch (last recipient of the gene if several gain events occurred).
	
	'refspetree' : 
	'ALEmodel'   : {'dated', 'undated'} specify the algorithm under which was generated the reconciliation (different parser)
	'verbose'    : integer (0-2)
	"""
	if method in ['last-gain', 'strict']:
		ogs, unclassified, dlabs = _prune_orthologs_bottom_up(recgt, ALEmodel=ALEmodel, **kw)
	elif method=='unicopy':
		ogs, dlabs = _prune_orthologs_top_down(recgt, ALEmodel=ALEmodel, **kw)
	elif method=='mixed':
		if 'gain_ogs' in kw:
			gain_ogs = kw['gain_ogs']
			dlabs = {}
		else:
			gain_ogs, unclassified, dlabs = _prune_orthologs_bottom_up(recgt, ALEmodel=ALEmodel, **kw)
		ogs, dlabs = _prune_orthologs_top_down(recgt, ALEmodel=ALEmodel, candidateOGs=gain_ogs, dlabs=dlabs, **kw)
	return ogs, dlabs

		
	
	
