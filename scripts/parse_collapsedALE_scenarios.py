#!/usr/bin/python

import os, sys, glob, getopt
import multiprocessing as mp
import cPickle as pickle
import shelve
import tree2
from ptg_utils import *
import parseALErec as pAr
import re

## Parameters
# file parsing parameter
constrainttag='mbconstraints'
replacedtag='-leaflabels_Spe2Pop.txt'

collapsedcladetag = '_CC-clade'
replacementcladetag = '_RC-clade'
replacementcladepat = re.compile('^.+_(.+_RC-clade[0-9]+)$')
outtaxlab = pAr.outtaxlab

# block reconstruction
gapsize = 2
gefagr = ['cds_code','replaced_cds_code','gene_family_id']

############ Functions

def eventLineages(recgt, dnodeallevt, ALEmodel='undated', recordEvTypes='DTS', onlyLeaves=[], deDupMatching=replacementcladepat):
	"""oreder events by gene lineage in the reconciled gene tree
	
	from a reconciled gene tree and the dict containing all the events in this tree in the format {node_id:(X, [don, ]rec), ...}, 
	return a list of event tuples for each gene lineage in the reconciled gene tree, i.e. the chain of events located above every tip of the tree"""
	def get_eventlineage(node, dnodeallevt, allevtlineages):
		if not node: return []
		nodeid = node.nodeid()
		eventpath = allevtlineages.get(nodeid)
		if not (eventpath is None):
			#~ print '#', nodeid
			#~ print 'found previous record'
			return eventpath
		else:
			# filter events
			eventpath = [evtup for evtup in dnodeallevt.get(nodeid, []) if (evtup[0] in recordEvTypes)]
			eventpath += get_eventlineage(node.father, dnodeallevt, allevtlineages)
			#~ print '#', nodeid
			#~ print 'allevtlineages:', allevtlineages
			#~ print 'eventpath:', eventpath
			if not node.is_leaf():
				#~ print 'recording'
				# record lineage of events from this node to the root ; dynamic programming !
				allevtlineages[nodeid] = eventpath
				#~ print 'allevtlineages:', allevtlineages
			return eventpath
	
	evtlineages = {}	# only lineages from the leaves to be returned
	allevtlineages = {} # cache dict for events at nodes shared by several leaves
	if ALEmodel=='undated':
		leavesandlabels = [(leaf, leaf.label().split('.')[0]) for leaf in recgt.get_leaves()]
	elif ALEmodel=='dated':
		#~ leavesandlabels = [(leaf, splitEventChain(leaf.label(), isleaf=True, ALEmodel='dated')[1]) for leaf in recgt.get_leaves()]
		leavesandlabels = [(leaf, leaf.label().split('.')[0].split('@')[0]) for leaf in recgt.get_leaves()]
	else:
		raise ValueError, "wrong ALE model specified: '%s'"%ALEmodel
	leavesandlabels.sort(key=lambda x: x[1]) # sort so that the representative first label to be captured by deDupMatching regex will be consistent across the recgt sample
	if onlyLeaves: leavesandlab = [x for x in leavesandlabels if (x[1] in onlyLeaves)]
	else: leavesandlab = leavesandlabels
	if deDupMatching:
		lela = []
		srctags = set([])
		for leaf, leaflab in leavesandlab:
			matchrc = deDupMatching.search(leaflab)
			if matchrc:
				rctag = matchrc.groups()[0]
				if not rctag in srctags:
					lela.append((leaf, leaflab))
					srctags.add(rctag)
			else:
				lela.append((leaf, leaflab))
		leavesandlab = lela
	for leaf, leaflab in leavesandlab:
		evtlineages[leaflab] = get_eventlineage(leaf, dnodeallevt, allevtlineages)
	return evtlineages

def translateRecStree(colspetree, refspetree):	
	"""matching branches of input collapsed species tree with those of full (i.e. uncollapsed) reference species tree; edit labels of collapsed tree and return dictionary of changed labels"""
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
	
def translateRecTs(nfrec, dcol2fullspenames, ALEmodel='undated'):
	"""edit input transfer report file, replacing collapsed species tree labels with full (i.e. uncollapsed) species treee labels; return dictionary of translated events and write out translated file"""
	devents = {}
	nfevent = nfrec.replace('ml_rec', 'Ts')
	fevent = open(nfevent, 'r')
	fevout = open(nfevent+'.fullScoords', 'w')
	for line in fevent:
		if ALEmodel=='undated':
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
	"""translate vector of events, replacing collapsed species tree labels with full (i.e. uncollapsed) species treee labels"""
	if len(ldtl)>0:
		if isinstance(str, ldtl[0]): return [dcol2fullspenames[loc] for loc in ldtl]
		else: return [tuple(dcol2fullspenames[x] for x in loc) for loc in ldtl]
	else:
		return ldtl
		
def translateEventLineage(deventlineages, dcol2fullspenames, drefspeeventTup2Ids=None, verbose=False):
	"""translate vectors of events stored in a dict by gene lineage, replacing collapsed species tree labels with full (i.e. uncollapsed) species treee labels
	
	if drefspeeventTup2Ids is provided, translated events will be stored using unique species tree event id (int) 
	instead of event tuple describing its type and location: (X, [don, ]rec).
	"""
	if dcol2fullspenames:
		trline = {}
		for nodelab, levtloc in deventlineages.iteritems():
			trltups = set()
			for evtloc in levtloc:
				trloc = tuple(dcol2fullspenames[x] for x in evtloc[1:])
				if len(trloc)>1 and len(set(trloc))==1:
					# trivial transfer event due to donor and the recipient in the collapsed species tree being nested in the full tree
					if verbose: print "ignore:", evtloc, ' ->', trloc
					continue
				trtup = evtloc[:1]+trloc
				trltups.add(trtup)
			trline[nodelab] = trltups
	else:
		trline = deventlineages
	if not drefspeeventTup2Ids:
		return trline
	else:
		# lighter version encoding events just by an integer referring to species tree events reference table
		return {nodelab:[drefspeeventTup2Ids[evtup] for evtup in levtup] for nodelab, levtup in trline.iteritems()}

def parseRec(nfrec, refspetree=None, ALEmodel='undated', drefspeeventTup2Ids=None, onlyLineages=[], recordEvTypes='DTS', minFreqReport=0, returnDict=True, lineageTableOutDir=None, noTranslateSpeTree=False, allEventByLineageByGenetree=False):
	"""parse reconciled gene tree sample, returning sampled events by gene lineage
	
	if allEventByLineageByGenetree is True, return more detailed data, stored in a dict with the following elements: 
	{
	 'allrectevtlineages': <dict of all single observed events by lineage by gene tree in the sample>, 
	 'devtlineagecount': <dict of all events and total observed frequency by lineage>, 
	 'dexactevt': <dict of frequencies of events, irrespective of the lineage in which they ocurred>'
	}
	otherwise (default), only the 'devtlineagecount' is returned.
	"""
	if not (returnDict or lineageTableOutDir): raise ValueError, "no output option chosen"
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
	if ALEmodel=='dated':
		# add reference for '#OUTSIDE#' taxon
		dcol2fullspenames[outtaxlab] = outtaxlab
	# parse reconciled gene trees
	# and extract (exact) event-wise event frequency
	dexactevt = {}
	devtlineagecount = {}
	allrectevtlineages = {}
	for i, recgt in enumerate(lrecgt):
		# gather scenario-scpecific events (i.e. dependent on reconciled gene tree topology, which varies among the sample)
		dlevt, dnodeallevt = pAr.parseRecGeneTree(recgt, colspetree, ALEmodel=ALEmodel, dexactevt=dexactevt, recgtsample=recgtsample, \
		                                          nsample=nsample, fillDTLSdict=False, recordEvTypes=recordEvTypes, \
		                                          excludeTaggedLeaves=collapsedcladetag, excludeTaggedSubtrees=replacementcladetag)
		# here events involving a replcement clade (RC) or leaf (CC) are excluded
		# * 'dexactevt' is used as cache to store frequencies of event s as inferred from regex searches of the event pattern
		# these frequencies are not specific to gene lineages, but aggregate the counts over the whole gene family
		# * 'dlevt' is of no use and here returned empty because of fillDTLSdict=False
		# would it not be empty, it could be translated to the full reference tree with:
		# tdlevt = {etype:translateEventList(ldtl, dcol2fullspenames, drefspeevents) for etype, ldtl in dlevt.iteritems()}
		evtlineages = eventLineages(recgt, dnodeallevt, ALEmodel=ALEmodel, onlyLeaves=onlyLineages, recordEvTypes=recordEvTypes)
		print 'evtlineages:', evtlineages
		tevtlineages = translateEventLineage(evtlineages, dcol2fullspenames, drefspeeventTup2Ids)
		print 'tevtlineages:', tevtlineages
		
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
			# but more efficient in memory use
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
		for geneleaflab, eventlineage in devtlineagecount.iteritems():
			for evtup, fevent in eventlineage.items():
				if float(fevent)/nsample < minFreqReport:
					del eventlineage[evtup]
	
	# optionally write out events gene by gene (those that occured at least once above a gene in [rooted] reconciled gene tree, and at which frequency)
	if lineageTableOutDir:
		nfTableEventsOut = os.path.join(lineageTableOutDir, "%s.%s.eventlineages"%(os.path.basename(nfrec), recordEvTypes))
		with open(nfTableEventsOut, 'w') as fTableOut:
			geneleaflabs = devtlineagecount.keys()
			geneleaflabs.sort()
			for geneleaflab in geneleaflabs:
				eventlineage = devtlineagecount[geneleaflab]
				for evtup, freq in eventlineage.iteritems():
					if drefspeeventTup2Ids: fTableOut.write('\t'.join((geneleaflab, str(evtup), str(freq)))+'\n')
					else: fTableOut.write('\t'.join((geneleaflab,)+evtup+(str(freq),))+'\n')
		print "stored events listed by gene lineage in '%s'"%nfTableEventsOut
	
	sys.stdout.flush()
	retd = {}
	retd['nfrec'] = nfrec
	if returnDict:
		retd['devtlineagecount'] = devtlineagecount
		if allEventByLineageByGenetree:
			retd['allrectevtlineages'] = allrectevtlineages
			retd['dexactevt'] = dexactevt
	else:
		return retd

# 
def parseRecTupArgs(args):
	"""wrapper function with arguments passed as a tuple"""
	nfrec, refspetree, ALEmodel, drefspeeventTup2Ids, onlyLineages, recordEvTypes, minFreqReport, returnDict, lineageTableOutDir = args
	return parseRec(nfrec, refspetree, ALEmodel, drefspeeventTup2Ids, onlyLineages, recordEvTypes, minFreqReport, returnDict, lineageTableOutDir)

def loadRecGeneTreeLabelAliasesAndListRecFiles(nflnfrec, nfgenefamlist=None, dircons=None, dirrepl=None, nbthreads=1, verbose=False):
	"""parse data relating to the genes and gene families to process.
	
	This includes the table of labels that were changed between original collapsed gene trees (as produced by script mark_unresolved_clades.py)
	and the species/population-tagged collapsed gene trees (as produced by script replace_species_by_pop_in_gene_trees.py)
	"""
	if nflnfrec:
		with open(nflnfrec, 'r') as flnfrec:
			lnfrec = [line.rstrip('\n') for line in flnfrec]
	else:
		lnfrec = []
	genefamlist = loadRecGeneTreeLabelAliases(nfgenefamlist, dircons=dircons, dirrepl=dirrepl, nbthreads=nbthreads, verbose=verbose)
	return (lnfrec, genefamlist)
	
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
	
def generateEventRefDB(refspetree, ALEmodel='undated', refTreeTableOutDir=None, TfromOutside=True):
	"""generates a dict of event tuples (of the form (X, [don, ]rec)). 
	
	Optionally writes out the reference species tree branches and event info to table files:
	- species_tree table dump has fields:       (branch_id, parent_branch_id, branch_name, is_tip)
	- species_tree_event table dump has fields: (event_id, event_type, don_branch_id, rec_branch_id)
	"""
	def recordGLEvent(eventid, evtup, outevtup, drefspeeventTup2Ids, drefspeeventId2Tups, foutspeevents=None):
		drefspeeventTup2Ids[evtup] = eventid
		drefspeeventId2Tups[eventid] = evtup
		if foutspeevents: foutspeevents.write('\t'.join([str(e) for e in outevtup])+'\n')
		return eventid+1
	
	drefspeeventTup2Ids = {}
	drefspeeventId2Tups = {}
	maxnid = 0
	eventid = 0
	if refTreeTableOutDir:
		foutspetree = open(os.path.join(refTreeTableOutDir, "phylogeny_species_tree.tab"), 'w')
		foutspeevents = open(os.path.join(refTreeTableOutDir, "phylogeny_species_tree_events.tab"), 'w')
	else:
		foutspeevents = None
	for node in refspetree:
		nid = node.nodeid()
		if nid>maxnid: maxnid = nid
		nlab = node.label()
		fnid = node.father_nodeid()
		#~ if ALEmodel=='undated':
		# for all models, do not record dates of events
		if refTreeTableOutDir: foutspetree.write('\t'.join((str(nid), str(fnid if fnid else ''), nlab, str(int(node.is_leaf()))))+'\n')
		#~ for et in pAr.eventTypes:
		for et in 'DTLS':
			evtup = (et, nlab)
			if et!='T':
				eventid = recordGLEvent(eventid, evtup, (eventid, et, '', nid), drefspeeventTup2Ids, drefspeeventId2Tups, foutspeevents)
			else:
				for donnode in refspetree:
					if not node is donnode:
						evtupt = evtup + (donnode.label(),)
						eventid = recordGLEvent(eventid, evtupt, (eventid, et, donnode.nodeid(), nid), drefspeeventTup2Ids, drefspeeventId2Tups, foutspeevents)
	
	# for backwards compatibility wwhen O event we not considered
	et = 'O'
	for node in refspetree:
		nid = node.nodeid()
		if nid>maxnid: maxnid = nid
		nlab = node.label()
		fnid = node.father_nodeid()
		evtup = (et, nlab)
		eventid = recordGLEvent(eventid, evtup, (eventid, et, '', nid), drefspeeventTup2Ids, drefspeeventId2Tups, foutspeevents)
	# add origination outside the tree
	outnid = maxnid + 1
	if refTreeTableOutDir: foutspetree.write('\t'.join((str(outnid), '', outtaxlab, '0'))+'\n')
	evtup = (et, outtaxlab)
	eventid = recordGLEvent(eventid, evtup, (eventid, et, '', outnid), drefspeeventTup2Ids, drefspeeventId2Tups, foutspeevents)
	
	if TfromOutside:
		et = 'T'
		# make this loop occur last so not to change the value of ids whether done or not
		for node in refspetree:
			nid = node.nodeid()
			nlab = node.label()
			evtupt = (et, outtaxlab, nlab)
			eventid = recordGLEvent(eventid, evtupt, (eventid, et, outnid, nid), drefspeeventTup2Ids, drefspeeventId2Tups, foutspeevents)
		
	if refTreeTableOutDir:
		foutspetree.close()
		foutspeevents.close()
	return (drefspeeventTup2Ids, drefspeeventId2Tups)

def parse_events(lnfrec, genefamlist=None, refspetree=None, ALEmodel='undated', \
                 drefspeeventTup2Ids={}, recordEvTypes='ODTS', minFreqReport=0, \
                 nfpickleEventsOut=None, nfshelveEventsOut=None, dirTableOut=None, nbthreads=1):
	"""from list of reconciliation files, families and genes to consider, return dictionary of reported events, by family and gene lineage"""
	lfams = [os.path.basename(nfrec).split('-')[0] for nfrec in lnfrec]
	if genefamlist: ingenes = [genefam.get('replaced_cds_code', genefam['cds_code']) for genefam in genefamlist if (genefam['gene_family_id'] in lfams)]
	else: ingenes=[]
	
	#~ # define wrapper function with set arguments, but the input reconciliation file
	#~ def parseRecSetArgs(nfrec):
		#~ return parseRec(nfrec, refspetree, drefspeeventTup2Ids, onlyLineages=ingenes, recordEvTypes=recordEvTypes, minFreqReport=minFreqReport, \
						#~ lineageTableOutDir=(os.path.join(dirTableOut, 'gene_tree_lineages') if dirTableOut else None))
	diroutab = (os.path.join(dirTableOut, 'gene_tree_lineages') if dirTableOut else None)
	returndict = bool(nfpickleEventsOut)
	iterargs = ((nfrec, refspetree, ALEmodel, drefspeeventTup2Ids, ingenes, recordEvTypes, minFreqReport, returndict, diroutab) for nfrec in lnfrec)
	
	# prepare output
	if nfshelveEventsOut:
		print "will gradually store event tuples in persistent dictionary (shelve) file: '%s'"%nfshelveEventsOut
		dfamevents = shelve.open(nfshelveEventsOut)
	else:
		dfamevents = {}
	
	if nbthreads==1:
		ildevents = (parseRecTupArgs(args) for args in iterargs)
	else:
		pool = mp.Pool(processes=nbthreads)
		#~ ldevents = pool.map(parseRecTupArgs, largs)
		# an iterator is returned by imap(); one needs to actually iterate over it to have the pool of parrallel workers to compute
		ildevents = pool.imap_unordered(parseRecTupArgs, iterargs)
	
	for deventfam in ildevents:
		fam = os.path.basename(deventfam['nfrec']).split('-')[0]
		devent = deventfam.get('devtlineagecount')
		dfamevents[fam] = devent
	
	if nfshelveEventsOut:
		print "saved 'dfamevents' to file '%s'"%nfpickleEventsOut
		dfamevents.close()
	elif nfpickleEventsOut:
		with open(nfpickleEventsOut, 'wb') as fpickleOut:
			pickle.dump(dfamevents, fpickleOut, protocol=2)
		print "saved 'dfamevents' to file '%s'"%nfpickleEventsOut
	
	return dfamevents

def main():

	opts, args = getopt.getopt(sys.argv[1:], 'hvT:', ['rec_sample_list=', 'ALE_algo=', \
	                                                'genefams=', 'dir_constraints=', 'dir_replaced=', \
	                                                'populations=', 'reftree=', \
	                                                'evtype=', 'minfreq=', \
	                                                'dir_table_out=', 'events_to_pickle=', 'events_to_shelve=', \
	                                                'threads=', 'help', 'verbose']) #, 'reuse=', 'max.recursion.limit=', 'logfile='
	dopt = dict(opts)
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	# reconciliation collection input
	nflnfrec = dopt['--rec_sample_list']
	
	# parsed events / matched events output options
	dirTableOut = dopt.get('--dir_table_out')
	nfpickleEventsOut = dopt.get('--events_to_pickle')
	nfshelveEventsOut = dopt.get('--events_to_shelve')
	if not (nfpickleEventsOut or nfshelveEventsOut or dirTableOut):
		raise ValueError, "an output option for parsed reconciliation must be chosen between '--dir_table_out', '--events_to_pickle' or '--events_to_shelve'"
	
	# other params
	# ALE reconciliation format
	ALEmodel = dopt.get('--ALE_algo', 'undated')
	
	# facultative input files
	nfpop = dopt.get('--populations')
	nfrefspetree = dopt.get('--reftree')
	nfgenefamlist = dopt.get('--genefams')
	dircons = dopt.get('--dir_constraints')
	dirrepl = dopt.get('--dir_replaced')
	
	# event filters
	recordEvTypes = dopt.get('--evtype', 'DTS')
	minFreqReport = float(dopt.get('--minfreq', 0))
	
	# runtime params
	nbthreads = int(dopt.get('--threads', dopt.get('-T', -1)))
	if nbthreads < 1: nbthreads = mp.cpu_count()
	verbose = ('-v' in dopt) or ('--verbose' in dopt)
	if verbose:
		print "# call: %s"%(' '.join(sys.argv))
		print "# dopt:", dopt
	
	if dirTableOut:
		ltd = ['ref_species_tree', 'gene_tree_lineages']
		for td in ltd:
			ptd = os.path.join(dirTableOut, td)
			if not os.path.isdir(ptd):
				os.mkdir(ptd)
	
	lnfrec, genefamlist = loadRecGeneTreeLabelAliasesAndListRecFiles(nflnfrec, nfgenefamlist, dircons, dirrepl, nbthreads=nbthreads)
	
	
	refspetree, dspe2pop = loadRefPopTree(nfrefspetree, nfpop)
	drefspeeventTup2Ids, drefspeeventId2Tups = generateEventRefDB(refspetree, ALEmodel, refTreeTableOutDir=(os.path.join(dirTableOut, 'ref_species_tree') if dirTableOut else None))
	
	dfamevents = parse_events(lnfrec, genefamlist, refspetree, ALEmodel, drefspeeventTup2Ids, recordEvTypes, minFreqReport, \
								  nfpickleEventsOut, nfshelveEventsOut, dirTableOut, nbthreads)


def usage():
	s = "Usage: [HELP MESSAGE INCOMPLETE]\n"
	s += "python %s --rec_sample_list /path/to/list_of_reconciliation_file_paths [OPTIONS]\n"%sys.argv[0]
	s += "Facultative options:\n"
	s += "\t\t--dir_constraints\tfolder containing files listing leaf labels of collapsed gene tree clades\n"
	s += "\t\t--dir_replaced\tfolder containing files listing replaced leaf labels (e.g. when giving a species identity to collapsed gene tree clades)\n"
	s += "\t\t--genefams\ttabulated file with header containing at least those two fields: 'cds_code', 'gene_family_id'\n"
	s += "\t\t\t\trows indicate the genes to be treated in the search, and to which gene family they belong\n"
	s += "\t\t\t\t(and hence in which reconciliation file to find them).\n"
	return s

################## Main execution

if __name__=='__main__':
	
	main()

