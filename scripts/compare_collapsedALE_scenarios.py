#!/usr/bin/python

import os, sys, glob, getopt
import multiprocessing as mp
import cPickle as pickle
import shelve
import itertools
import gc
import numpy as np
from ptg_utils import *
# try the power of Cython
import pyximport
pyximport.install()
#~ pyximport.install(pyimport=True) # cythonize all the python modules avalaible
from coevol_score import coevol_lineages

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

def _select_lineage_event_clause_factory(evtypes, valtoken, lineagetable):
	rlocdsIJ = "INNER JOIN %s USING (replacement_label_or_cds_code)"%lineagetable
	if evtypes:
		evtyperestrictIJ = "INNER JOIN species_tree_events USING (event_id)"
		evtyperestrictWC = "AND event_type IN %s"%repr(tuple(e for e in evtypes))
	else:
		evtyperestrictWC = evtyperestrict = ""
	return (rlocdsIJ, evtyperestrictIJ, evtyperestrictWC)

def _select_lineage_event_query_factory(aBYb, evtypes, valtoken, lineagetable, \
                                        joinTable=None, addselcols=(), distinct=False, \
                                        operator='=', addWhereClause='', orderBy=''):
	a, b = aBYb
	rlocdsIJ, evtyperestrictIJ, evtyperestrictWC = _select_lineage_event_clause_factory(evtypes, valtoken, lineagetable)
	if joinTable:
		bIJ = 'INNER JOIN %s USING (%s)'%(joinTable, b)
		bWC = ''
	else:
		bIJ = ''
		bWC = 'AND %s%s%s'%(b, operator, valtoken)
	ob = "ORDER BY %s"%orderBy if orderBy else ''
	w = 'WHERE' if (evtyperestrictWC or addWhereClause or bWC) else ''
	tqfields = (('DISTINCT' if distinct else ''), repr((a,)+tuple(addselcols)).strip('(,)').replace("'", ''), 
				rlocdsIJ, evtyperestrictIJ, bIJ,
				w, evtyperestrictWC, bWC, addWhereClause, ob)
	preq = "SELECT %s %s FROM gene_lineage_events %s %s %s %s %s %s %s %s ;"%tqfields
	return preq.replace('WHERE AND ', 'WHERE ')

def _query_create_temp_events_lineage(gene, preq, dbcur, temptablename):
	creq = "CREATE TEMP TABLE %s AS "%temptablename + preq
	dbcur.execute(creq, (gene,)) 

def _query_events_lineage(gene, preq, dbcur, with_create_temp=None):
	if with_create_temp:
		creq = "CREATE TEMP TABLE %s AS "%with_create_temp + preq
		dbcur.execute(creq, (gene,)) 
		dbcur.execute("SELECT * FROM %s ;"%with_create_temp)
	else:
		dbcur.execute(preq, (gene,)) 
	return dbcur.fetchall()

def _query_events_lineage_sorted(gene, preq, dbcur, sortfield=0, with_create_temp=None):
	return sorted(_query_events_lineage(gene, preq, dbcur, with_create_temp), key=lambda x: x[sortfield])

def _query_matching_lineage_event_profiles(args, timing=False, verbose=False):
	
	if timing: time = __import__('time')
	querylineageidfam, dbname, dbengine, nsample, evtypes, baseWC, diffamWC, minevjointfreq, lineagetable = args
	querylineage_id, queryfam = querylineageidfam
	if verbose: print 'querylineage_id:', querylineage_id, queryfam
	dbcon, dbcul, dbtype, valtoken = get_dbconnection(dbname, dbengine)
	# first get the vector of (event_id, freq) tuples in focal lineage
	preq_evbyli = _select_lineage_event_query_factory(('event_id', 'rlocds_id'), \
	                                                  evtypes, valtoken, lineagetable, \
	                                                  addselcols=('freq as f0',), \
	                                                  addWhereClause=baseWC)
	if verbose: print preq_evbyli%querylineage_id
	tempeventtable = 'events_of_rlocds_id%d'%querylineage_id
	#~ lineage_eventfreqs = _query_events_lineage_sorted(querylineage_id, preq_evbyli, dbcul, with_create_temp=tempeventtable)
	_query_create_temp_events_lineage(querylineage_id, preq_evbyli, dbcul, tempeventtable)
	#~ dlineage_eventfreqs = dict(lineage_eventfreqs)
	#~ lineage_events = tuple(ef[0] for ef in lineage_eventfreqs)
	#~ if not lineage_events: return []
	#~ if verbose: print lineage_events
	if dbcul.rowcount == 0:
		# no events associated withthis lineage, return empty array
		# -- this sort of pointless query should be avoided by
		# not listing for query lineages without selectable events
		# i.e. pre-filter input to this function call
		return []
	
	# then build a list of gene lineages to compare, based on common occurence of at least N event (here only 1 common event required)
	# and filtering by {same|different|all} gene families
	famWC = diffamWC%queryfam if '%' in diffamWC else diffamWC
	# and the id of compared lineage to be > reference lineage, to avoid duplicate comparisons
	# or equal (for same-family scope mode), i.e. self-comparison, to evaluate the the maximum association score for this lineage
	lineageorderWC=" AND rlocds_id >= %d"%querylineage_id
	preq_libyev = _select_lineage_event_query_factory(('rlocds_id', 'event_id'), \
	                                                  evtypes, valtoken, lineagetable, \
	                                                  addselcols=('f0', 'freq as f1',), \
	                                                  joinTable=tempeventtable, \
	                                                  addWhereClause=famWC+baseWC+lineageorderWC, orderBy='rlocds_id')
	
	if verbose: print preq_libyev
	dbcul.execute(preq_libyev) 
	nsamplesq = nsample**2
	coevollineages = coevol_lineages(dbcul, querylineage_id, nsamplesq, minevjointfreq)
	
	dbcul.execute("drop table %s ;"%tempeventtable)
	dbcon.close()
	if verbose: print coevollineages
	return coevollineages

def dbquery_matching_lineage_event_profiles(dbname, dbengine='postgres', \
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
		dbcur.execute("""create table %s as select DISTINCT replacement_label_or_cds_code, gene_family_id, rlocds_id 
		                  from replacement_label_or_cds_code2gene_families 
		                  inner join ingenes using (replacement_label_or_cds_code) ;
		              """%(lineagetable, ))
		dbcur.execute("drop table ingenes;")
		dbcon.commit()
	
	# prepare invariant segments of by-lineage queries
	minevfWC = " AND gene_lineage_events.freq >= %d"%int(mineventfreq*nsample) if mineventfreq>0 else ''
	maxevfWC = " AND gene_lineage_events.freq  < %d"%int(maxeventfreq*nsample) if maxeventfreq<1 else ''
	baseWC = minevfWC+maxevfWC
	if matchScope=='between_fams':
		diffamWC = " AND gene_family_id != '%s'" # %fam to plugged in in child function call
	elif matchScope=='within_fams':
		diffamWC = " AND gene_family_id = '%s'" # %fam to plugged in in child function call
	elif matchScope=='all':
		diffamWC = ""
	else:
		raise ValueError, "incorrect value '%s' for variable 'matchScope'"%repr(matchScope)
	def make_arg_tup(lineage_id):
		#~ return (lineage_id, dbname, dbengine, nsample, evtypes, matchScope, mineventfreq, maxeventfreq, minevjointfreq, lineagetable)
		return (lineage_id, dbname, dbengine, nsample, evtypes, baseWC, diffamWC, minevjointfreq, lineagetable)
		
	# get set of gene lineages
	qlibyev = "SELECT rlocds_id, gene_family_id FROM %s;"%lineagetable
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
		#~ matchesOutRad = os.path.basename(nfgenefamlist).rsplit('.', 1)[0] if nfgenefamlist else ''
		matchesOutRad = ''
		matchesOutDirRad = os.path.join(dirTableOut, 'gene_lineage_assocations', matchScope+'_scores', matchesOutRad)
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
