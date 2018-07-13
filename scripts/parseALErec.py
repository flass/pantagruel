#!/usr/bin/python
import tree2
import sys, os, getopt
import re

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
		if node.is_leaf() and (excludeTaggedLeaves in nodelab):
			# the species assignment of this leaf is not certain (inferred)
			# and thus events leading directly to this leaf are not to be trusted and should be skipped
			return
		nodeid = node.nodeid()
		if not nodelab:
			print node
			raise ValueError, "unannotated node"
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
def findNextDatedEventSlice(nodelab, receptionOK=False):
	nextev1 = nodelab.find('@')
	if nextev1 > 0:
		if nodelab[nextev1-1] in {'T', 'D'}:
			nextev1 -=1		# account for the preceding D or T
		elif not receptionOK:
			raise IndexError
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
	"""yet to implement"""
	raise NotImplementedError
	lineage = nodelab.split('.') # bring in code from parseUndatedRecGeneTree()
	return lineage

def getLeafLabelDatedEvent(nodelab):
	# remove leaf label prefix before parsing (possible) trailing event chain
	nextev = findNextDatedEventSlice(nodelab, receptionOK=True)
	leaflab = nodelab[:nextev] # the whole string if nextev is None
	return leaflab

def splitSingleDatedEvent(nodelab, isleaf=False):
	"""from a gene tree branch event chain, yield the type and location of every event, and if a leaf, finishes by the leaf label"""
	if isleaf:
		# remove leaf label prefix before parsing (possible) trailing event chain
		nextev = findNextDatedEventSlice(nodelab, receptionOK=True)
		leaflab = nodelab[:nextev] # the whole string if nextev is None
		# fist yield
		yield leaflab
		s = nodelab[len(leaflab):]
	else:
		s = nodelab
	i = 0
	while s:
		print s
		reception = False
		# process separators in order
		if s.startswith('.T@'):
			# transfer-loss // emission of transfer by donor followed by loss of resident copy; just changes species identity of the gene tree branch (from donor to recipent species)
			# will be followed by transfer reception event
			reception = True
			thisev = 3
			evtype = 'TdL'
		elif s.startswith('.'):
			if i==0:
				# simple speciation
				evtype = 'S'
			else:
				# speciation-loss // speciation followed by loss in the sister lineage; just changes species identity of the gene tree branch (from parent to daughter species)
				evtype = 'SL'
			thisev = 1
		elif s.startswith('@'):
			# reception of transfer by recipient = gain
			thisev = 1
			evtype = 'Tr'
		elif s.startswith('T@'):
			# emission of transfer by donor; not relevant to gene content
			thisev = 2
			evtype = 'Td'
		elif s.startswith('D@'):
			# duplication = gain
			thisev = 2
			evtype = 'D'
		elif isleaf:
			# final speciation
			evtype = 'S'
			thisev = 0
		else:
			raise ValueError
		i += 1
		nextev = findNextDatedEventSlice(s[thisev:], receptionOK=reception)
		if nextev is None:
			# no other event
			yield (evtype, s[thisev:])
			s = None
		else:
			yield (evtype, s[thisev:nextev+thisev])
			s = s[nextev+thisev:]

def splitEventChain(nodelab, isleaf=False, ALEmodel='dated'):
	if ALEmodel=='undated':
		lineage = [event for event in splitSingleDatedEvent(nodelab, isleaf=isleaf)]
	elif ALEmodel=='dated':
		lineage = [event for event in splitSingleDatedEvent(nodelab, isleaf=isleaf)]
		if isleaf:
			# remove the first-yielded leaf label from vent list
			leaflab = lineage[0]
			lineage = lineage[1:]
		else:
			leaflab = None
		# line of events is read left-to-right FORWARD in time; reverse it for a BACKWARD result
		lineage = reversed(lineage)
	return (lineage, leaflab)

def _prune_orthologs_bottom_up(node, ALEmodel='dated', **kw):
	"""strict definition of ortholous groups: only those related by a line of speciation events (no transfer, no duplication) are orthologs"""
	# initiate at root
	unclassified = set([])
	orthologGroups = kw.get('orthologGroups', [])
	dlabs = kw.get('dlabs', {})
	nodelab = node.label()
	if not nodelab: raise ValueError, "unannotated node:\n%s"%str(node)
	lineage, leaflab = splitEventChain(nodelab, isleaf=node.is_leaf(), ALEmodel=ALEmodel)
	if node.is_leaf():
		unclassified.add(leaflab) # where unclasified set is filled up
		dlabs[nodelab] = leaflab
	else:
		# post-order traversal recursion:
		for child in node.children:
			# first explore children
			orthologGroups, uncl, dlabs = _prune_orthologs_bottom_up(child, orthologGroups=orthologGroups, dlabs=dlabs, ALEmodel=ALEmodel)
			unclassified |= uncl
	# then explore node itself
	if not unclassified:
		# all leaves below have been already classified in orthologous group; stop here.
		return orthologGroups, unclassified, dlabs, orthoSubs
	else:
		# list of events must go BACKWARD in time when read left-to-right
		for event in lineage:
			evtype, evloc = event
			if evtype in {'Tr', 'T'}:
				# gene was last gained here by a species ancestor:
				# what was not yet classifed in this subtree is an orthologous group
				ortho = tuple(sorted(unclassified))
				orthologGroups.append(ortho)
				unclassified = set([])
				break
			elif evtype=='D':
				# gene was duplicated here by species ancestor:
				# what was not yet classifed in the child subtrees
				# of this subtree are each an orthologous group
				for child in node:
					subuncl = set([dlabs[leaflabevchain] for leaflabevchain in child.iter_leaf_labels()]) & unclassified
					ortho = tuple(sorted(subuncl))
					orthologGroups.append(ortho)
				unclassified = set([])
				break
		return orthologGroups, unclassified, dlabs

def _prune_orthologs_top_down(node, refspetree=None, ALEmodel='dated', **kw):
	assert refspetree
	# initiate at root
	orthologGroups = kw.get('orthologGroups', [])
	dlabs = kw.get('dlabs', {})
	splitparam0, splitparam1 = kw.get('splitparam', ('_',1))
	# first establish the dictionary of actual leaf lables, without the trailing event chain
	for nodelab in node.iter_leaf_labels():
		dlabs[nodelab] = getLeafLabelDatedEvent(nodelab)
	if not nodelab: raise ValueError, "unannotated node:\n%s"%str(node)
	lineage, leaflab = splitEventChain(nodelab, isleaf=node.is_leaf(), ALEmodel=ALEmodel)
	# list of events must go BACKWARD in time when read left-to-right
	for event in lineage:
		evtype, evloc = event
		if evtype in {'Tr', 'T'}:
			leaflabs = [dlabs[leaflabevchain] for leaflabevchain in node.iter_leaf_labels()]
			lspe = [leaflab.split(splitparam0, splitparam1)[0] for leaflab in leaflabs]
			sspe = set(lspe)
			if len(lspe) < len(sspe):
				# not unicopy
				break # for event loop
			else:
				speanc = refspetree[evloc]
				assert isinstance(speanc, tree2.Node)
				clademrca = refspetree.mrca(sspe)
				if clademrca.is_childorself(speanc):
					ortho = tuple(sorted(leaflabs))
					orthologGroups.append(ortho)
					break # for event loop
	else:
		# if did not find the clade to be an orthologous group, recurse down the tree:
		for child in node.children:
			orthologGroups, dlabs = _prune_orthologs_top_down(node, orthologGroups=orthologGroups, dlabs=dlabs, refspetree=refspetree, ALEmodel=ALEmodel)
	return orthologGroups, dlabs

def getOrthologues(recgt, ALEmodel='dated', method='strict', **kw):
	"""prune reconciled gene tree at transfer or duplication branches to identify subs of orthologues. 
	
	several methods are possible:
	'strict':
	    bottom-up approach (tip-to-root) where only those leaves related by a line of speciation events (no transfer, no duplication) are orthologs; 
	    the search is achieved by 
	'unicopy_clades':
	    top-down approach (root-to-tips) where orthologous groups are defined as the gene tree clades containing only unique representatives of species  
	    that belong to the species clade below the ancestor to which a gene gain is mapped (last recipient of the gene if several gain events occurred 
	    on the gene tree branch). This allows to recognize as orthologous groups  clades where transfer events within the clade that do not induce 
	    change of copy number (when considering the sum of all events in the subtree), i.e. gene conversion events
	'mixed': (heuristic ?) approach mixing the two previous: pruning of transfers in the bottom-up approach may reveal bona-fide unicopy clades of orthologs.
	
	using the last two options requires to provide a reference species tree, passed on with the keyword argument 'refspetree'.
	"""
	if method=='strict':
		return _prune_orthologs_bottom_up(recgt, ALEmodel=ALEmodel, **kw)
	elif method=='unicopy_clades':
		refspetree = kw['refspetree']
		return _prune_orthologs_top_down(recgt, refspetree=None, ALEmodel=ALEmodel, **kw)
	


		
	
	
