#!/usr/bin/python

import os, sys, glob, getopt
import multiprocessing as mp
import cPickle as pickle
import itertools
import numpy as np
import tree2
from replace_species_by_pop_in_gene_trees import annotatePopulationInSpeciesTree, parseMrBayesConstraints, getdspe2pop
from mark_unresolved_clades import mean, median
import parseALErec as pAr


############ Functions

def eventLineages(recgt, dnodeallevt, onlyLeaves=[]):
	def get_eventlineage(node, dnodeallevt, allevtlineages):
		if not node: return []
		nodeid = node.nodeid()
		if nodeid in allevtlineages: return allevtlineages[nodeid]
		# record lineage of events from this node to the root ; dynamic programming !
		eventpath = dnodeallevt.get(nodeid, []) + get_eventlineage(node.father, dnodeallevt, allevtlineages)
		allevtlineages[nodeid] = eventpath
		return eventpath
	
	evtlineages = {}	# only lineages from the leaves to be returned
	allevtlineages = {} # cache dict for events at nodes shared by several leaves
	leavesandlabels = [(leaf, leaf.label().split('.')[0]) for leaf in recgt.get_leaves()]
	if onlyLeaves: leavesandlab = [x for x in leavesandlabels if (x[1] in onlyLeaves)]
	else: leavesandlab = leavesandlabels
	for leaf, leaflab in leavesandlab:
		evtlineages[leaflab] = get_eventlineage(leaf, dnodeallevt, allevtlineages)
	return evtlineages

def translateRecStree(colspetree, refspetree):	
	"""edit labels of input collapsed species tree in reference to fully labelled, uncollapsed species tree and return dictionary of changed labels"""
	dcol2fullspenames = {}
	tcolspetree = colspetree.deepcopybelow()
	for node in tcolspetree:
		if node.is_leaf():
			dcol2fullspenames[node.label()] = node.label()
		else:
			coresplab = refspetree.coalesce(node.get_leaf_labels()).label()
			dcol2fullspenames[node.label()] = coresplab
			node.edit_label(coresplab)
	return [tcolspetree, dcol2fullspenames]
	
def translateRecTs(nfrec, dcol2fullspenames, mode='undated'):
	devents = {}
	nfevent = nfrec.replace('ml_rec', 'Ts')
	fevent = open(nfevent, 'r')
	fevout = open(nfevent+'.fullScoords', 'w')
	for line in fevent:
		if mode=='undated':
			if line.startswith('#'):
				fevout.write(line)
				continue
			lsp = line.rstrip('\n').split('\t')
			#~ print lsp
			tdon = dcol2fullspenames[lsp[1]]
			trec = dcol2fullspenames[lsp[2]]
			fevout.write('\t'.join([tdon, trec, lsp[3]])+'\n')
			devents.setdefault(trec, {})[tdon] = float(lsp[3])
	fevent.close()
	fevout.close()
	#~ print devents
	return devents

def translateEventList(ldtl, dcol2fullspenames):
	if len(ldtl)>0:
		if isinstance(str, ldtl[0]): return [dcol2fullspenames[loc] for loc in ldtl]
		else: return [tuple(dcol2fullspenames[x] for x in loc) for loc in ldtl]
	else:
		return ldtl
		
def translateEventLineage(deventlineages, dcol2fullspenames, drefspeeventTup2Ids=None, verbose=False):
	trline = {}
	for nodelab, levtloc in deventlineages.iteritems():
		trltups = []
		for evtloc in levtloc:
			trloc = tuple(dcol2fullspenames[x] for x in evtloc[1:])
			if len(trloc)>1 and len(set(trloc))==1:
				# trivial transfer event due to donor and the recipient in the collapsed species tree being nested in the full tree
				if verbose: print "ignore:", evtloc, ' ->', trloc
				continue
			trtup = evtloc[:1]+trloc
			trltups.append(trtup)
		trline[nodelab] = trltups
	if not drefspeeventTup2Ids:
		return trline
	else:
		# lighter version encoding events just by an integer referring to species tree events reference table
		return {nodelab:[drefspeeventTup2Ids[evtup] for evtup in levtup] for nodelab, levtup in trline.iteritems()}

def parseRec(nfrec, refspetree=None, drefspeeventTup2Ids=None, onlyLineages=[], recordEvTypes='T', minFreqReport=0, lineageTableOutDir=None, noTranslateSpeTree=False, allEventByLineageByGenetree=False):
	print nfrec
	# parse reconciliation file and extract collapsed species tree, mapping of events (with freq.) on the species tree, and reconciled gene trees
	colspetree, subspetree, lrecgt, recgtlines, restrictlabs, dnodeevt = pAr.parseALERecFile(nfrec)
	nsample = len(lrecgt)
	recgtsample = ''.join(recgtlines)
	if not noTranslateSpeTree:
		tcolspetree, dcol2fullspenames = translateRecStree(colspetree, refspetree)
	else:
		# no need to translate
		tcolspetree, dcol2fullspenames = colspetree, {}
		if refspetree:
			assert refspetree.hasSameTopology(tcolspetree, checkInternalLabels=True)
	# parse reconciled gene trees
	# and extract (exact) event-wise event frequency
	dexactevt = {}
	devtlineagecount = {}
	allrectevtlineages = {}
	for i, recgt in enumerate(lrecgt):
		# gather scenario-scpecific events (i.e. dependent on reconciled gene tree topology, which varies among the sample)
		dlevt, dnodeallevt = pAr.parseRecGeneTree(recgt, colspetree, dexactevt, recgtsample, nsample, fillDTLSdict=False, recordEvTypes=recordEvTypes)
		# * 'dexactevt' is used as cache to store frequencies of event s as inferred from regex searches of the event pattern
		# these frequencies are not specific to gene lineages, but aggregate the counts over the whole gene family
		# * 'dlevt' is of no use and here returned empty because of fillDTLSdict=False
		# would it not be empty, it could be translated to the full reference tree with:
		# tdlevt = {etype:translateEventList(ldtl, dcol2fullspenames, drefspeevents) for etype, ldtl in dlevt.iteritems()}
		evtlineages = eventLineages(recgt, dnodeallevt, onlyLeaves=onlyLineages)
		tevtlineages = translateEventLineage(evtlineages, dcol2fullspenames, drefspeeventTup2Ids)
		
		if allEventByLineageByGenetree:
			# one way to proceed is to build the object 'allrectevtlineages'
			# a dict that contains all events in a lineage, 
			# for all the lineages in reconciled gene tree, 
			# for all the reconcile gene trees in the ALE sample.
			# IT CAN BE A VERY HEAVY OBJECT.
			for geneleaflab, evtlineage in tevtlineages.iteritems():
				allrectevtlineages.setdefault(geneleaflab, []).append(evtlineage)
		else:
			# another way is to aggregate data immediately
			# might be slower due to many updates of the 'devtlineagecount' dict,
			# but mor efficient in memory use
			for geneleaflab, evtlineage in tevtlineages.iteritems():
				for evtup in evtlineage:
					nevtup = devtlineagecount.setdefault(geneleaflab, {}).setdefault(evtup, 0)
					devtlineagecount[geneleaflab][evtup] = nevtup + 1
	
	if allEventByLineageByGenetree:
		devtlineagecount = {}
		for geneleaflab, allreclineages in allrectevtlineages.iteritems():
			allrecevt = reduce(lambda x, y: x+y, allreclineages)
			# combine event counts across the sample
			fevent = {evtup:allrecevt.count(evtup) for evtup in set(allrecevt)}
			if minFreqReport>0:
				# skips the low-frequency events
				if float(fevent)/nsample < minFreqReport: continue
			devtlineagecount[geneleaflab] = fevent
	elif minFreqReport>0:
		# cleanup by deleting low-frequency events a posteriori
		for geneleaflab in devtlineagecount:
			for evtup, fevent in devtlineagecount[geneleaflab].iteritems():
				if float(fevent)/nsample < minFreqReport:
					del devtlineagecount[geneleaflab][evtup]
	
	# optionally write out events gene by gene (those that occured at least once above a gene in [rooted] reconciled gene tree, and at which frequency)
	if lineageTableOutDir:
		nfTableOut = os.path.join(lineageTableOutDir, "%s.%s.eventlineages"%(os.path.basename(nfrec), recordEvTypes))
		with open(nfTableOut, 'w') as fTableOut:
			geneleaflabs = devtlineagecount.keys()
			geneleaflabs.sort()
			for geneleaflab in geneleaflabs:
				eventlineage = devtlineagecount[geneleaflab]
				for evtup, freq in eventlineage.iteritems():
					if drefspeeventTup2Ids: fTableOut.write('\t'.join((geneleaflab, str(evtup), str(freq)))+'\n')
					else: fTableOut.write('\t'.join((geneleaflab,)+evtup+(str(freq),))+'\n')
		print "stored events listed by gene lineage in '%s'"%nfTableOut
	
	if allEventByLineageByGenetree:
		retd = {}
		retd['allrectevtlineages'] = allrectevtlineages
		retd['devtlineagecount'] = devtlineagecount
		retd['dexactevt'] = dexactevt
		return retd
	else:
		return devtlineagecount
			
def loadNamesAndListRecFiles(nflnfrec, nfgenefamlist, dircons, dirrepl, verbose=False):
	with open(nflnfrec, 'r') as flnfrec:
		lnfrec = [line.rstrip('\n') for line in flnfrec]
	with open(nfgenefamlist, 'r') as fgenefamlist:
		header = fgenefamlist.readline().strip(' \n').split(',')
		# remove trailing whitespace present after gene_family_id ; account for PostgreSQL footer line count of fetched rows
		genefamlist = [dict(zip(header, line.replace(' ,', ',').strip(' \n').split(','))) for line in fgenefamlist if (not line.endswith('rows)\n'))] 
	lfams = [os.path.basename(nfrec).split('-')[0] for nfrec in lnfrec]
	dreplacedlab = {}
	for fam in lfams:
		globcons = "%s/%s/*%s*"%(dircons, fam, constrainttag)
		if verbose: print globcons
		lnfcons = glob.glob(globcons)
		constraintclades = parseMrBayesConstraints(lnfcons)
		globreplaced = "%s/%s*%s"%(dirrepl, fam, replacedtag)
		if verbose: print globreplaced
		lnfreplaced = glob.glob(globreplaced)
		if len(lnfreplaced)!=1:
			sys.stderr.write("Warning! will use first file among those found matching pattern '%s'\nas reference for collapsed gene tree tip label aliases:\n%s\n"%( globreplaced, '\n'.join(lnfreplaced) ))
		nfreplaced = lnfreplaced[0]
		with open(nfreplaced, 'r') as freplaced:
			for line in freplaced:
				claorgene, repllab = line.rstrip('\n').split('\t')
				if claorgene in constraintclades:
					for genelab in constraintclades[claorgene]:
						dreplacedlab[genelab] = repllab
				else:
					dreplacedlab[claorgene] = repllab
	
	for genefam in genefamlist:
		genelab = genefam['cds_code']
		genrep = dreplacedlab.get(genelab, genelab)
		genefam['replaced_cds_code'] = genrep
	return [lnfrec, lfams, genefamlist]
	
def loadRefPopTree(nfrefspetree, nfpop):
	# annotate reference species tree with ancestral population node labels
	lnamepops = []
	with open(nfpop, 'r') as fpop:
		for line in fpop:
			if line.startswith('#'): continue
			lsp = line.rstrip('\n').split('\t')
			lnamepops.append((lsp[0], tuple(lsp[1].split())))
	refspetree = tree2.AnnotatedNode(file=nfrefspetree)
	refspetree.complete_internal_labels(order=0, ffel=True)
	refspetree.complete_node_ids(force=True)
	annotatePopulationInSpeciesTree(refspetree, lnamepops, returnCopy=False, returnAncNodes=False)
	dspe2pop = getdspe2pop(lnamepops)
	nfrefspetreeout = nfrefspetree.rsplit('.', 1)[0]+'_internalPopulations.nwk'
	refspetree.write_newick(nfrefspetreeout, ignoreBS=True)
	return (refspetree, dspe2pop)
	
def generateEventRefDB(refspetree, model='undated', refTreeTableOutDir=None):
	"""generates a dict of event tuples (of the form (X, rec, [don, ])). 
	
	Optionally writes out the reference species tree branches and event info to table files
	"""
	drefspeeventTup2Ids = {}
	drefspeeventId2Tups = {}
	eventid = 0
	if refTreeTableOutDir:
		foutspetree = open(os.path.join(refTreeTableOutDir, "phylogeny_species_tree.tab"), 'w')
		foutspeevents = open(os.path.join(refTreeTableOutDir, "phylogeny_species_tree_events.tab"), 'w')
	for node in refspetree:
		nid = node.nodeid()
		nlab = node.label()
		fnid = node.father_nodeid()
		if model=='undated':
			if refTreeTableOutDir: foutspetree.write('\t'.join((str(nid), str(fnid if fnid else ''), nlab, str(int(node.is_leaf()))))+'\n')
			for et in pAr.eventTypes:
				evtup = (et, nlab)
				if et!='T':
					drefspeeventTup2Ids[evtup] = eventid
					drefspeeventId2Tups[eventid] = evtup
					if refTreeTableOutDir: foutspeevents.write('\t'.join((str(eventid), et, str(nid), ''))+'\n')
					eventid += 1
				else:
					for donnode in refspetree:
						if not node is donnode:
							evtupt = evtup + (donnode.label(),)
							drefspeeventTup2Ids[evtupt] = eventid
							drefspeeventId2Tups[eventid] = evtupt
							if refTreeTableOutDir: foutspeevents.write('\t'.join([str(e) for e in (eventid, et, nid, donnode.nodeid())])+'\n')
							eventid += 1
	if refTreeTableOutDir:
		foutspetree.close()
		foutspeevents.close()
	return (drefspeeventTup2Ids, drefspeeventId2Tups)

def matchEventInLineages(dfamevents, genei, fami, genej, famj, noSameFam=True, blocks={}, dspe2pop={}, drefspeeventId2Tups={}, eventtypes='T', excludeNodeLabels=[], **kw):
	"""Generator function yielding matched homologous events from different reconciled scenarios. Proceeds by dissecting gene tree tip-to-root lineages.
	
	It typically involes comparing scenarios from different gene families (enforced if noSameFam=True), 
	but could aslo in theory be ued to compare scenarios from the same gene family, 
	under different reconciliation inference settings (e.g. different DTL rates, 
	or different reference tree collapsing). 
	
	These families can have different reference species tree, or more specifically, 
	different collapsed version of the same full reference species tree.
	
	This function takes a input a dict-of-dict-of dict object 'dfamevents', containing the 
	lists all the sampled event lineages, sorted by families and genes: 
	dfamevents[fami][genei] should thus return a dict of the following form: 
	{(X, rec, [don, ]):freq, ...}, 
	where keys are tuples of DTL event types ('X' in {D,T,L,S}) and addresses, defined as 'rec' and optionally 'don' 
	(for horizontal transfers), the labels in species tree of node where the event occurred,
	OR, if drefspeeventId2Tups is provided:
	{i:freq, ...}
	where i is the species tree event unique reference id,
	and values are 'freq' the event frequency within the sample of N reconciled gene trees.
	
	!!! Because closely related genes share lineage history that are duplicated in this representation,
	!!! this is a highly redundant representation and thus quite heavy format.
	
	This function is a generator and thus yields matched event tuples 
	
	If 'dspe2pop' is provided, alias addresses in the compared refence species tree will be generated,
	so that events involving strains/lineage from within a population (in a species tree where it is not collapsed) 
	and those involving the collapsed population can be recognised as homologous.
	
	However, it is though safer to not try and match events occurring below the population ancestor nodes,
	because of the low reolution of the population subtree making inferences to be incertain, 
	and to a lesser extant to match events occurring below the population ancestor nodes, 
	due to the possibly poor quality of the heuristic used to tag populations on collapsed gene tree clades.
	
	Such events of uncertain qulity can be excluded from the comparison using 'excludeNodeLabels' argument.
	"""
	def eventLocSpeOrPop(evtup, dspe2pop):
		return [evtup[:1]+tuple(y) for y in itertools.product(*[list(set([ x, dspe2pop.get(x, x) ])) for x in evtup[1:]])]
	
	if noSameFam and fami==famj:
		# stop immediately
		return
	
	devi = dfamevents[fami][genei] #; print devi
	devj = dfamevents[famj][genej] #; print devj
	klause = genei=='ESCCOL13_5175' and genej=='ECOS51_ENTCGC014269_CC002'
	evtupis = kw.get('events', devi.keys())
	#~ for evtupi in evtupis:
		#~ # every event in genei lineage should be match at most once in genej lineage
		#~ if evtupi[0] not in eventtypes: continue
		#~ for akaevtupi in eventLocSpeOrPop(evtupi, dspe2pop):
			#~ for evtupj in devj:
				#~ akaevtupjs = eventLocSpeOrPop(evtupj, dspe2pop)
				#~ if akaevtupi in akaevtupjs:
					#~ if genej in blocks.get(evtupi, []): continue
					#~ yield [evtupi, devi[evtupi], devj[evtupj]]
					#~ break # the for evtupj loop
	for evtupi in evtupis:
		# one event tuple per sampled recongiled gene tree, representing the gene lineage (= tip-to-root node path in gene tree)
		if klause : print 'evtupi:', evtupi
		# every event in genei lineage should be matching at most once in genej lineage
		matched = False
		if drefspeeventId2Tups.get(evtupi, evtupi)[0] not in eventtypes: continue
		if dspe2pop: akaevtupis = eventLocSpeOrPop(evtupi, dspe2pop)
		else: akaevtupis = (evtupi,)
		if klause : print 'akaevtupis:', akaevtupis 
		for evtupj in devj:
			if klause:  print 'evtupj:', evtupj 
			if matched: break # the for evtupj loop
			# mapping of aliases: species to collapsed population
			if dspe2pop: akaevtupjs = eventLocSpeOrPop(evtupj, dspe2pop)
			else: akaevtupjs = (evtupj,)
			if klause:  print 'akaevtupjs:', akaevtupjs 
			for akaevtupi in akaevtupis:
				# one event lineage tuple per sampled recongiled gene tree
				if akaevtupi in akaevtupjs:
					if klause : print 'nosino'
					if not (genej in blocks.get(evtupi, [])):
						# avoid yielding the event match again if was previously spotted in a gene block
						if klause : print 'rhhhoooooo', [evtupi, devi[evtupi], devj[evtupj]]
						yield [evtupi, devi[evtupi], devj[evtupj]]
					matched = True
					# only one alias should be matched (as with alias vs. alias matching, it can give up to 2^len(evloc) matches)
					break # the for akaevtupi loop
					
	#~ if klause: raise ValueError
	
	
def searchPWMatches(dfamevents, genefams, nfoutblocks, nfpicklePWMatches, dspe2pop={}, drefspeeventId2Tups={}, eventtypes='T'):
	foutblocks = open(nfoutblocks, 'w')
	foutblocks.write('\t'.join(['event', 'don', 'rec', 'len.block', 'sum.freq', 'min.freq', 'max.freq', 'mean.freq', 'median.freq', 'genes', 'eventFreqs'])+'\n' )
	PWMatches = {}
	for i in range(len(genefams)-1):
		genei, genri, fami = [genefams[i][f] for f in gefagr]
		print genei, genri, fami
		for j in range(i+1, len(genefams)):
			genej, genrj, famj = [genefams[j][f] for f in gefagr]
			#~ print genej, genrj, famj
			for matchedT, fi, fj in matchEventInLineages(dfamevents, genri, fami, genrj, famj, dspe2pop=dspe2pop, drefspeeventId2Tups=drefspeeventId2Tups, eventtypes=eventtypes):
				PWMatches.setdefault((genei, genej), []).append((matchedT, fi, fj)) 
				
	with open(nfpicklePWMatches, 'wb') as fpicklePWMatches:
		 pickle.dump(PWMatches, fpicklePWMatches, protocol=2)
	print "saved 'PWMatches' to file '%s'"%nfpicklePWMatches
	return PWMatches				

def treeExtantEventSet(levents, refspetree, fun=max):
	"""measure dispersion of a set of events on the species tree
	
	by default return maximum distance, can be average (fun=mean), minimum (fun=min)...
	"""
	ldist = []
	if len(levents) <=1: return maxdist
	# from the set of all (recipient) locations, compute the maximum distance on the species tree
	# ... should do that on the consensus gene tree... if one could map each sampled reconciled gene tree to the consensus!!!
	allTreeLocs = list(set(reduce(lambda x,y: x[-1:]+y[-1:], levents)))
	for loci in allTreeLocs:
		ni = refspetree[loci]
		for locj in allTreeLocs:
			if loci!=locj:
				nj = refspetree[locj]
				d = ni.distance(nj)
				ldist.append(d)
	if fun is max: ldist.append(0)
	return fun(ldist)

def reconstructBlocks(dfamevents, genefams, nfoutblocks, nfpicklePWMatches, dspe2pop={}, gapsize=-1):
	"""finds blocks of consecutive genes that share homologous events"""
	foutblocks = open(nfoutblocks, 'w')
	foutblocks.write('\t'.join(['event', 'don', 'rec', 'len.block', 'sum.freq', 'min.freq', 'max.freq', 'mean.freq', 'median.freq', 'genes', 'eventFreqs'])+'\n' )
	dblockTs = {}
	dblockTsReplGenes = {}
	dblockTfreqs = {}
	PWMatches = {}
	for i in range(len(genefams)-1):
		j = i+1
		genei, genri, fami = [genefams[i][f] for f in gefagr]
		genej, genrj, famj = [genefams[j][f] for f in gefagr]
		#print "%s: %s (%s) <-> %s: %s (%s)"%(fami, genei, genri, famj, genej, genrj)
		for matchedT, fi, fj in matchEventInLineages(dfamevents, genri, fami, genrj, famj, dspe2pop=dspe2pop, blocks=dblockTsReplGenes):
			#print matchedT, 
			dblockTs.setdefault(matchedT, [genei]).append(genej)
			dblockTsReplGenes.setdefault(matchedT, [genri]).append(genrj)
			dblockTfreqs.setdefault(matchedT, [fi]).append(fj)
			#PWMatches.setdefault((genei, genej), []).append((matchedT, fi, fj)) 
			print genei, genej, matchedT, fi, fj
			# try to extend transfer block
			k = j+1
			moreTs = True
			genel = genej ; fl = fj
			if gapsize>=0: curgap = gapsize
			else: curgap = 1
			while (k<len(genefams) and (moreTs or curgap>0)):
				genek, genrk, famk = [genefams[k][f] for f in gefagr]
				moreTs = [Tff for Tff in matchEventInLineages(dfamevents, genri, fami, genrk, famk, dspe2pop=dspe2pop, blocks=dblockTsReplGenes, events=[matchedT])]
				if moreTs:
					curgap = gapsize
					fk = moreTs[0][2]
					#print genek,
					dblockTs[matchedT].append(genek)
					dblockTsReplGenes[matchedT].append(genrk)
					dblockTfreqs[matchedT].append(fk)
					#PWMatches.setdefault((genel, genek), []).append((matchedT, fl, fk)) 
					print ' '*k, genel, genek, matchedT, fl, fk
					genel = genek ; fl = fk
				else:
					if gapsize>=0: curgap -= 1
				k += 1
			else:
				print ''
				if gapsize>=0:
					# remove occurrence of theat transfer pattern for further independent creation of another block
					blockTgenes = dblockTs.pop(matchedT)
					blockTgenesR = dblockTsReplGenes.pop(matchedT)
					blockTfreqs = dblockTfreqs.pop(matchedT)
				else:
					# gap of infinite size ; just test co-transfer linkage without considering gene neighborhood
					blockTgenes = dblockTs[matchedT]
					blockTgenesR = dblockTsReplGenes[matchedT]
					blockTfreqs = dblockTfreqs.pop[matchedT]
				for x, gx in enumerate(blockTgenes):
					for y, gy in enumerate(blockTgenes):
						if x==y: continue
						PWMatches.setdefault((gx, gy), []).append((matchedT, blockTfreqs[x], blockTfreqs[y])) 
				
				foutblocks.write('\t'.join(list(matchedT) \
				+ [str(fun(blockTfreqs)) for fun in (len, sum, min, max, mean, median)] \
				+ [' '.join(blockTgenes), ' '.join([str(f) for f in blockTfreqs])])+'\n' )

	foutblocks.close()
	with open(nfpicklePWMatches, 'wb') as fpicklePWMatches:
		 pickle.dump(PWMatches, fpicklePWMatches, protocol=2)
	print "saved 'PWMatches' to file '%s'"%nfpicklePWMatches
	return PWMatches

def usage():
	s = ""
	s += "\t\t--dir_replaced\tfolder containing files listing replaced leaf labels (e.g. wen collapsing gene tree clades)\n"
	s += "\t\t--genefams\ttabulated file with header containing at least those two fields: 'cds_code', 'gene_family_id'\n"
	s += "\t\t\t\trows indicated the genes to be treated in the search, and to which gene family they belong (and hence in which reconciliation file to find them)\n"
	return s

################## Main execution

## Parameters
# file parsing parameter
constrainttag='mbconstraints'
replacedtag='-leaflabels_Spe2Pop.txt'

# block reconstruction
gapsize = 2
gefagr = ['cds_code','replaced_cds_code','gene_family_id']

# stat summary param
threshMinEventFreq = 10
PWMstats = [ \
'PWMatchesSummedFreq']#, \
#~ 'PWMatchesUniqueCount', \
#~ 'PWMatchesUniqueCountMinEvFq', \
#~ 'PWMatchesTreeExtent', \
#~ 'PWMatchesTreeExtentMinEvFq']

if __name__=='__main__':

	opts, args = getopt.getopt(sys.argv[1:], 'hv', ['rec_sample_list=', 'pickled_events=', 'dir_constraints=', 'dir_replaced=', \
	                                               'populations=', 'reftree=', 'genefams=', \
	                                               'evtype=', 'minfreq=', 'dir_table_out=', \
												   'threads=', 'help', 'verbose']) #, 'reuse=', 'verbose=', 'max.recursion.limit=', 'logfile='
	dopt = dict(opts)
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	try:
		nflnfrec = dopt['--rec_sample_list']
		loaded = False
	except KeyError:
		try:
			nfpickle = dopt['--pickled_events']
			with open(nfpickle, 'rb') as fpickle:
				dfamevents = pickle.load(fpickle)
			print "loaded 'dfamevents' from file '%s'"%nfin
			nflnfrec = nfpickle.split(pickletag)[0]
			loaded = True
		except KeyError:
			raise ValueError, "must provide input file through either '--rec_sample_list' or '--pickled_events' options"
	
	dircons = dopt['--dir_constraints']
	dirrepl = dopt['--dir_replaced']
	nfpop = dopt['--populations']
	nfrefspetree = dopt['--reftree']
	nfgenefamlist = dopt['--genefams']
	recordEvTypes = dopt.get('--evtype', 'DTS')
	minFreqReport = float(dopt.get('--minfreq', 0))
	dirTableOut = dopt.get('--dir_table_out',)
	nbthreads = int(dopt.get('--threads', 1))
	verbose = ('-v' in dopt) or ('--verbose' in dopt)
	
	pickletag = '.parsedRecs.%s.pickle'%recordEvTypes
	
	if dirTableOut:
		for td in ['ref_species_tree', 'gene_tree_lineages']:
			ptd = os.path.join(dirTableOut, td)
			if not os.path.isdir(ptd):
				os.mkdir(ptd)
	
	refspetree, dspe2pop = loadRefPopTree(nfrefspetree, nfpop)
	drefspeeventTup2Ids, drefspeeventId2Tups = generateEventRefDB(refspetree, model='undated', refTreeTableOutDir=(os.path.join(dirTableOut, 'ref_species_tree') if dirTableOut else None))
	
	lnfrec, lfams, genefamlist = loadNamesAndListRecFiles(nflnfrec, nfgenefamlist, dircons, dirrepl)

	if not loaded:	
		ingenes = [genefam['replaced_cds_code'] for genefam in genefamlist if (genefam['gene_family_id'] in lfams)]
		# define wrapper function with set arguments, but the input reconciliation file
		def parseRecSetArgs(nfrec):
			return parseRec(nfrec, refspetree, drefspeeventTup2Ids, onlyLineages=ingenes, recordEvTypes=recordEvTypes, minFreqReport=minFreqReport, \
			                lineageTableOutDir=(os.path.join(dirTableOut, 'gene_tree_lineages') if dirTableOut else None))

		if nbthreads==1:
			ldevents = [parseRecSetArgs(nfrec) for nfrec in lnfrec]
		else:
			pool = mp.Pool(processes=nbthreads)
			ldevents = pool.map(parseRecSetArgs, lnfrec)
			
		dfamevents = dict(zip(lfams, ldevents))
		
		nfpickle = nflnfrec+pickletag
		with open(nfpickle, 'wb') as fpickle:
			pickle.dump(dfamevents, fpickle, protocol=2)
		print "saved 'dfamevents' to file '%s'"%nfpickle
	##
	genefams = [genefam for genefam in genefamlist if (genefam['gene_family_id'] in lfams)]

	nfpicklePWMatches = nfgenefamlist.rsplit('.', 1)[0]+'.PWgeneEventMatches.%s.pickle'%recordEvTypes
	if os.access(nfpicklePWMatches, os.F_OK):
		with open(nfpicklePWMatches, 'rb') as fpicklePWMatches:
			PWMatches = pickle.load(fpicklePWMatches)
		print "loaded 'PWMatches' from file '%s'"%nfpicklePWMatches
	else:
		PWMatches = searchPWMatches(dfamevents, genefams, nfoutblocks, nfpicklePWMatches, dspe2pop=dspe2pop, drefspeeventId2Tups=drefspeeventId2Tups, eventtypes=recordEvTypes)
	# summarize
	genelist = [genefam['cds_code'] for genefam in genefams]

	PWMmats = {PWMstat:np.zeros((len(genelist), len(genelist)), dtype=float) for PWMstat in PWMstats}
	for genepair, levff in PWMatches.iteritems():
		genepairi = tuple(genelist.index(geni) for geni in genepair)
		pwmatches = [ev for ev, fi, fj in levff]
		PWMmats['PWMatchesSummedFreq'][genepairi] = float(sum(reduce(lambda x,y: x[1:]+y[1:], levff, (0,))))/2
		#~ pwmatchesMinEvFq = [ev for ev, fi, fj in levff if (fi>=threshMinEventFreq or fj>=threshMinEventFreq)]
		#~ PWMmats['PWMatchesUniqueCount'][genepairi] = len(pwmatches) 
		#~ PWMmats['PWMatchesUniqueCountMinEvFq'][genepairi] = len(pwmatchesMinEvFq) 
		#~ PWMmats['PWMatchesTreeExtent'][genepairi] = treeExtantEventSet(pwmatches, refspetree)
		#~ PWMmats['PWMatchesTreeExtentMinEvFq'][genepairi] = treeExtantEventSet(pwmatchesMinEvFq, refspetree)

	# write output
	for PWMstat, PWMmat in PWMmats.iteritems():
		nfmat = nfgenefamlist.rsplit('.', 1)[0]+'.%s.%s.csv'%(PWMstat, recordEvTypes)
		with open(nfmat, 'w') as fmat:
			fmat.write(','+','.join(genelist)+'\n')
			#~ PWMmat.tofile(fmat, sep=',')
			for i, gene in enumerate(genelist):
				fmat.write(','.join([gene]+[str(x) for x in PWMmat[i,]])+'\n')
			print "wrote %s output into file '%s'"%(PWMstat, nfmat)
		
