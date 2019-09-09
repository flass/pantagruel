#!/usr/bin/python

import os, sys, glob, getopt
from parse_collapsedALE_scenarios import *

def parseRec(nfrec, refspetree=None, recformat='ecceTERA', drefspeeventTup2Ids=None, onlyLineages=[], recordEvTypes='DTS', minFreqReport=0, returnDict=True, \
             lineageTableOutDir=None, noTranslateSpeTree=False, allEventByLineageByGenetree=False, verbose=False):
	"""parse reconciled gene tree sample, returning sampled events by gene lineage
	
	if allEventByLineageByGenetree is True, return more detailed data, stored in a dict with the following elements: 
	{
	 'allrectevtlineages': <dict of all single observed events by lineage by gene tree in the sample>, 
	 'devtlineagecount': <dict of all events and total observed frequency by lineage>, 
	 'dexactevt': <dict of frequencies of events, irrespective of the lineage in which they ocurred>'
	}
	otherwise (default), only the 'devtlineagecount' is returned.
	"""
	pass

def main():

	opts, args = getopt.getopt(sys.argv[1:], 'hvT:', ['rec_sample_list=', 'rec_format=', \
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
	# ecceTERA reconciliation format: either native ecceTERA format (.txt extension) or Mowgli format (.mr extension)
	recformat = dopt.get('--rec_format', 'ecceTERA')
	
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
								  nfpickleEventsOut, nfshelveEventsOut, dirTableOut, nbthreads, verbose)


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

