#!/usr/bin/python

import os, sys, glob, getopt
import multiprocessing as mp
import cPickle as pickle
import shelve
import itertools
import gc
import numpy as np
from ptg_utils import *

## Parameters
# block reconstruction
gapsize = 2
gefagr = ['cds_code','replaced_cds_code','gene_family_id']

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

# here there would be scope for Cython static defs as this is done a lot
def mulBoundInt2Float(a, b, scale):
	f1 = float(a)/scale
	f2 = float(b)/scale
	return f1*f2

def _select_lineage_event_clause_factory(evtypes, valtoken, lineagetable):
	rlocdsIJ = "INNER JOIN %s USING (replacement_label_or_cds_code)"%lineagetable
	if evtypes:
		evtyperestrictIJ = "INNER JOIN species_tree_events USING (event_id)"
		evtyperestrictWC = "AND event_type IN %s"%repr(tuple(e for e in evtypes))
	else:
		evtyperestrictWC = evtyperestrict = ""
	return (rlocdsIJ, evtyperestrictIJ, evtyperestrictWC)

def _select_lineage_event_query_factory(aBYb, evtypes, valtoken, lineagetable, joinTable=None, addselcols=(), distinct=False, operator='=', addWhereClause=''):
	a, b = aBYb
	rlocdsIJ, evtyperestrictIJ, evtyperestrictWC = _select_lineage_event_clause_factory(evtypes, valtoken, lineagetable)
	if joinTable:
		bIJ = 'INNER JOIN %s USING (%s)'%(joinTable, b)
		bWC = ''
	else:
		bIJ = ''
		bWC = '%s%s%s'%(b, operator, valtoken)
	w = 'WHERE' if (evtyperestrictWC or addWhereClause or bWC) else ''
	tqfields = (('DISTINCT' if distinct else ''), repr((a,)+tuple(addselcols)).strip('(,)').replace("'", ''), 
				rlocdsIJ, evtyperestrictIJ, bIJ,
				w, evtyperestrictWC, bWC, addWhereClause)
	preq = "SELECT %s %s FROM gene_lineage_events %s %s %s %s %s %s %s ;"%tqfields
	return preq.replace('WHERE AND ', 'WHERE ')

def _query_events_by_lineage(gene, preq, dbcur, with_create_temp=None):
	if with_create_temp:
		creq = "create temp table %s as "%with_create_temp + preq
		dbcur.execute(creq, (gene,)) 
		dbcur.execute("select * from %s ;"%with_create_temp)
	else:
		dbcur.execute(preq, (gene,)) 
	return dbcur.fetchall()

def _query_events_by_lineage_sorted(gene, preq, dbcur, sortfield=0, with_create_temp=None):
	return sorted(_query_events_by_lineage(gene, preq, dbcur, with_create_temp), key=lambda x: x[sortfield])

def _query_matching_lineage_event_profiles(args, use_gene_labels=False, timing=False, arraysize=10000, verbose=False):
	
	if timing: time = __import__('time')
	lineageidfam, dbname, dbengine, nsample, evtypes, baseWC, minevjointfreq, lineagetable = args
	lineage_id, fam = lineageidfam
	if verbose: print 'lineage_id:', lineage_id, fam
	dbcon, dbcul, dbtype, valtoken = get_dbconnection(dbname, dbengine)
	# first get the vector of (event_id, freq) tuples in lineage
	preq_evbyli = _select_lineage_event_query_factory(('event_id', 'rlocds_id'), \
	                                                  evtypes, valtoken, lineagetable, \
	                                                  addselcols=('freq',), \
	                                                  addWhereClause=minevfWC+maxevfWC)
	if verbose: print preq_evbyli%lineage_id
	templineagetable = 'lineages_rlocds_id%d'%lineage_id
	lineage_eventfreqs = _query_events_by_lineage_sorted(lineage_id, preq_evbyli, dbcul, with_create_temp=templineagetable)
	dlineage_eventfreqs = dict(lineage_eventfreqs)
	lineage_events = tuple(ef[0] for ef in lineage_eventfreqs)
	if not lineage_events: return []
	if verbose: print lineage_events
	
	# then build a list of gene lineages to compare, based on common occurence of at least N event (here only 1 common event required)
	# and filtering by {same|different|all} gene families
	basefamWC = baseWC%fam if '%' in baseWC else baseWC
	# and the id of compared lineage to be > reference lineage, to avoid duplicate comparisons
	lineageorderWC=" AND rlocds_id > %d"%lineage_id
	preq_libyev = _select_lineage_event_query_factory(('rlocds_id', 'event_id'), \
	                                                  evtypes, valtoken, lineagetable, \
	                                                  joinTable=templineagetable, distinct=True, \
	                                                  addWhereClause=basefamWC+lineageorderWC)
	if verbose: print preq_libyev
	dbcul.execute(preq_libyev) 
	match_lineages = dbcul.fetchall()
	dbcul.execute("drop table %s ;"%templineagetable)
	
	lmatches = []
	for tmatch_lineage_id in match_lineages:
		jointfreq = 0.0
		match_lineage_id = tmatch_lineage_id[0]
		mg_lineage_eventfreqs = _query_events_by_lineage(match_lineage_id, preq_evbyli, dbcul)
		for eid, f in mg_lineage_eventfreqs:
			f0 = dlineage_eventfreqs.get(eid)
			#~ if f0: jointfreq += mulBoundInt2Float(f0, f, nsample, maxval=1.0)
			if f0: jointfreq += mulBoundInt2Float(f0, f, nsample)
		if jointfreq >= minevjointfreq:
			lmatches.append( (lineage_id, match_lineage_id, jointfreq) )
	
	dbcon.close()
	if verbose: print lmatches
	return lmatches

def dbquery_matching_lineage_event_profiles(dbname, dbengine='postgres', use_gene_labels=False, \
                                            genefamlist=None, exclRecSpeBranches=[], matchScope='between_fams', \
                                            nsample=1.0, evtypes=None, mineventfreq=0.0, maxeventfreq=1.0, minevjointfreq=0.0, \
                                            matchesOutDirRad=None, nfpickleMatchesOut=None, returnList=False, \
                                            nbthreads=1, verbose=False, **kw):
	
	def output_match_line(lm, lmatches, fout, nfoutrad, kfout, foutMaxSize=1024**3):
		# check if output file max size has been reached
		if fout: 
			foutsize = fout.tell()
			if foutsize >= foutMaxSize:
				fout.close()
				kfout +=1
				fout = open(nfoutrad+'.%d'%kfout, 'w')
		if (nfpickleMatchesOut or returnList): lmatches += lm
		if matchesOutDirRad or verbose:
			for tgpcf in lm:
				matchline = '%d\t%d\t%f\n'%tgpcf
				if fout: fout.write(matchline)
			if verbose:
				print len(lm), 'matches'
				sys.stdout.flush()
		return (fout, kfout)
	
	timing = kw.get('timing')
	if timing: time = __import__('time')							
	dbcon, dbcur, dbtype, valtoken = get_dbconnection(dbname, dbengine)
	nsamplesq = nsample*nsample
	if not (matchesOutDirRad or nfpickleMatchesOut or returnList):
		raise ValueError, "must specify at least one output option among: 'matchesOutDirRad', 'nfpickleMatchesOut', 'returnList'"
	
	lineagetable = 'replacement_label_or_cds_code2gene_families'
	if genefamlist:
		lineagetable += '_restricted'
		# create real table so can be seen by other processes
		ingenes = []
		for genefam in genefamlist:
			genetreelab = genefam.get('replacement_label_or_cds_code', genefam.get('replaced_cds_code', genefam.get('cds_code')))
			if genetreelab: ingenes.append( (genetreelab,) )
			else: raise KeyError, "no appropriate field defining the gene tree label in %s"%repr(genefam)
		dbcur.execute("drop table if exists ingenes;")
		dbcur.execute("create table ingenes (replacement_label_or_cds_code VARCHAR(60));" )
		dbcur.executemany("insert into ingenes values (%s) ;"%(valtoken), ingenes)
		dbcur.execute("drop table if exists %s;"%lineagetable)
		dbcur.execute("""create table %s as select replacement_label_or_cds_code, gene_family_id, rlocds_id 
		                  from replacement_label_or_cds_code2gene_families 
		                  inner join ingenes using (replacement_label_or_cds_code) ;
		              """%(lineagetable, ))
		dbcur.execute("drop table ingenes;")
		dbcon.commit()
	
	# prepare invariant segments of by-lineage queries
	minevfWC = " AND gene_lineage_events.freq >= %d"%int(mineventfreq*nsample) if mineventfreq>0 else ''
	maxevfWC = " AND gene_lineage_events.freq  < %d"%int(maxeventfreq*nsample) if maxeventfreq<1 else ''
	if matchScope=='between_fams':
		diffamWC = " AND gene_family_id != '%s'" # %fam to plugged in in child function call
	elif matchScope=='within_fams':
		diffamWC = " AND gene_family_id = '%s'" # %fam to plugged in in child function call
	elif matchScope=='all':
		diffamWC = ""
	else:
		raise ValueError, "incorrect value '%s' for variable 'matchScope'"%repr(matchScope)
	baseWC = diffamWC+minevfWC+maxevfWC
	def make_arg_tup(lineage_id):
		#~ return (lineage_id, dbname, dbengine, nsample, evtypes, matchScope, mineventfreq, maxeventfreq, minevjointfreq, lineagetable)
		return (lineage_id, dbname, dbengine, nsample, evtypes, baseWC, minevjointfreq, lineagetable)
		
	# get set of gene lineages
	qlibyev = "select rlocds_id, gene_family_id from %s;"%lineagetable
	if verbose: print qlibyev
	dbcur.execute(qlibyev)
	ltlineageidfams = dbcur.fetchall()
				
	
	if matchesOutDirRad:
		nfoutrad = os.path.join(matchesOutDirRad, 'matching_events.tab')
		kfout = 0
		fout = open(nfoutrad+'.%d'%kfout, 'w')
	else:
		kfout = nfoutrad = fout = None
	lmatches = []
	if nbthreads==1:
		for i, tlineageidfam in enumerate(ltlineageidfams):
			lm = _query_matching_lineage_event_profiles(make_arg_tup(tlineageidfam), verbose=max(verbose-1, 0))
			fout, kfout = output_match_line(lm, lmatches, fout, nfoutrad, kfout)
			if verbose: sys.stdout.write("\r%d\t"%i)
	else:
		pool = mp.Pool(processes=nbthreads)
		iterargs = (make_arg_tup(tlineageidfam) for tlineageidfam in ltlineageidfams)
		iterlm = pool.imap_unordered(_query_matching_lineage_event_profiles, iterargs, chunksize=1)
		# an iterator is returned by imap_unordered(); one needs to actually iterate over it to have the pool of parrallel workers to compute
		for lm in iterlm:
			fout, kfout = output_match_line(lm, lmatches, fout, nfoutrad, kfout)
	
	if nfpickleMatchesOut:
		with open(nfpickleMatchesOut, 'wb') as fpickleOut:
			pickle.dump(lmatches, fpickleOut, protocol=pickle.HIGHEST_PROTOCOL)
			# can be a more efficient option for disk space than a simple table,
			# but the required writing time and the accumulated memory space might make it redibitory
		print "saved 'lmatches' to file '%s'"%nfpickleMatchesOut
	if matchesOutDirRad:
		fout.close()
	dbcon.close()
	if returnList:
		return lmatches

def main():
	# options relating to:
	#  - scenario format
	#  - event matching filters and scope
	#  - event dataset input
	#  - matched events output
	#  - program runtime
	opts, args = getopt.getopt(sys.argv[1:], 'hv', ['ALE_algo=', 'nrec_per_sample=', \
	                                                'exclude_species_tree_branches=', 'event_type=', 'min_freq=', 'max_freq=', 'min_joint_freq=', 'match_scope=', \
	                                                'genefams=', 'dir_constraints=', 'dir_replaced=', \
	                                                'events_from_pickle=', 'events_from_shelve=', 'events_from_postgresql_db=', 'events_from_sqlite_db=', \
	                                                'matches_to_shelve=', 'dir_table_out=', \
	                                                'threads=', 'help', 'verbose=']) #, 'reuse=', 'max.recursion.limit=', 'logfile='
	dopt = dict(opts)
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	matchScope = dopt.get('--match_scope', 'between_fams')
	goodscopes = ['between_fams', 'within_fams', 'all']
	if matchScope not in goodscopes:
		raise ValueError, "valid values for --match_scope argument are: %s"%repr(goodscopes).strip('[]')
	
	# reconciliation collection / parsed events input options
	nfpickleEventsIn = dopt.get('--events_from_pickle')
	nfshelveEventsIn = dopt.get('--events_from_shelve')
	dbname   = dopt.get('--events_from_postgresql_db', dopt.get('--events_from_sqlite_db'))
	if nfshelveEventsIn:
		loaddfamevents = 'shelve'
	elif nfpickleEventsIn:
		loaddfamevents = 'pickle'
	elif dbname and ('--events_from_postgresql_db' in dopt):
		loaddfamevents = 'postgres'
	elif dbname and ('--events_from_sqlite_db' in dopt):
		loaddfamevents = 'sqlite'
	else:
		raise ValueError, "must provide input file through either '--events_from_pickle', '--events_from_shelve', '--events_from_sqlite_db' or '--events_from_postgresql_db' options"
	
	# parsed events / matched events output options
	dirTableOut = dopt.get('--dir_table_out')
	nfpickleMatchesOut = dopt.get('--matches_to_pickle')
	nfshelveMatchesOut = dopt.get('--matches_to_shelve')
	
	# other params
	
	# normalization factor (and max per lineage) of observed event frequencies
	nrecsample = dopt.get('--nrec_per_sample', 1000.0)
	# facultative input files
	nfgenefamlist = dopt.get('--genefams')
	dircons = dopt.get('--dir_constraints')
	dirrepl = dopt.get('--dir_replaced')
	
	# event filters
	recordEvTypes = dopt.get('--event_type', 'DTS')
	minFreqReport = float(dopt.get('--min_freq', 0.0))
	maxFreqReport = float(dopt.get('--max_freq', 1.0))
	minJointFreqReport = float(dopt.get('--min_joint_freq', 0.0))
	exclRecSpeBranches = str(dopt.get('--exclude_species_tree_branches', '')).split(',')
	# runtime params
	nbthreads = int(dopt.get('--threads', 1))
	verbose = int(dopt.get('--verbose', dopt.get('-v', 0)))
	
	if dirTableOut:
		ltd = ['gene_event_matches']
		for td in ltd:
			ptd = os.path.join(dirTableOut, td)
			if not os.path.isdir(ptd):
				os.mkdir(ptd)
	
	genefamlist = loadRecGeneTreeLabelAliases(nfgenefamlist, dircons, dirrepl, nbthreads=nbthreads, verbose=verbose)		
	
	if dirTableOut:
		matchesOutRad = os.path.basename(nfgenefamlist).rsplit('.', 1)[0] if nfgenefamlist else ''
		matchesOutDirRad = os.path.join(dirTableOut, 'gene_event_matches', matchesOutRad)
		if not os.path.isdir(matchesOutDirRad):
			os.mkdir(matchesOutDirRad)
	else:
		matchesOutDirRad = None
	dbquery_matching_lineage_event_profiles(dbname, dbengine=loaddfamevents, genefamlist=genefamlist, nsample=nrecsample, \
                         evtypes=recordEvTypes, exclRecSpeBranches=exclRecSpeBranches, matchScope=matchScope, \
                         mineventfreq=minFreqReport, maxeventfreq=maxFreqReport, minevjointfreq=minJointFreqReport, \
                         nfpickleMatchesOut=nfpickleMatchesOut, nfshelveMatchesOut=nfshelveMatchesOut, matchesOutDirRad=matchesOutDirRad, \
                         nbthreads=nbthreads, verbose=verbose)

def usage():
	s = "Usage: [HELP MESSAGE INCOMPLETE]\n"
	s += "python %s --events_from_{pickle|shelve|postgresql_db|sqlite_db} source [--dir_table_out dest] [--matches_to_{pickle|shelve} dest] [OTHER OPTIONS]\n"%sys.argv[0]
	s += "Facultative options:\n"
	s += "\t\t--dir_constraints\tfolder containing files listing leaf labels of collapsed gene tree clades\n"
	s += "\t\t--dir_replaced\tfolder containing files listing replaced leaf labels (e.g. when giving a species identity to collapsed gene tree clades)\n"
	s += "\t\t--genefams\ttabulated file with header containing at least those two fields: 'cds_code', 'gene_family_id'\n"
	s += "\t\t\t\trows indicate the genes to be treated in the search, and to which gene family they belong\n"
	s += "\t\t\t\t(and hence in which reconciliation file to find them).\n"
	s += "\t\t--match_scope={'between_fams'|'within_fams'|'all'} define scope where co-evolution score is to be evaluated:\n"
	s += "\t\t\t\t- within or between gene families only, or both. Between families is the default behavioue;\n"
	s += "\t\t\t\t- within_families can be chosen to allow to cluster closely related lineages with significantly shared ancestry within families\n"
	s += "\t\t\t\t  and to restrict accordingly the search for matches between families\n."
	return s

################## Main execution

if __name__=='__main__':
	
	main()
