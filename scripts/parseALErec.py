#!/usr/bin/python
import tree2
import sys, os, getopt
import re

eventTypes = 'DTLS'

def getOriSpeciesFromEventLab(eventlab, sgsep='_'):
	# split at DT location separator '@', then possibly at T don/rec seprator '->', and finally shorten the species label if node is a leaf
	elab = eventlab.split('@')[1] if '@' in eventlab else eventlab
	return elab.split('->')[0].split(sgsep)[0]

def parseRecGeneTree(recgt, spet, dexactevt, recgtsample, nsample, sgsep='_', restrictlabs=[], fillDTLSdict=True, recordEvTypes='DTL', excludeTaggedLeaves=None, excludeTaggedSubtrees=None):
	"""extract list of events from one reconciled gene tree as found in output of ALEml_undated (Szolosi et al., 2013; https://github.com/ssolo/ALE)
	
	This function returns two objects:
	
	- a dict object 'dlevt', containing the (possibly redundant) list of unique events having occurred in the tree,
	sorted by event type, which can be 'D', 'T', 'L' or 'S', for duplications, horizontal transfers, losses or speciations.
	This record is a crude aggregate of what happened in the whole gene family, present for back-compatiility purposes.
	The behaviour of filling up this record can be turned off to gain time and memory with fillDTLSdict=False.
	
	- a dict object 'dnodeallevt', containing the list all events in one reconciliation scenario
	(a single scenario anomg the sample), sorted by gene tree node id. 
	Event are represented as tuples of the following form: (X, rec, [don, ] freq),
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
			alllabext = [lab.split('_', 1)[1] for lab in node.get_leaf_labels()]
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
