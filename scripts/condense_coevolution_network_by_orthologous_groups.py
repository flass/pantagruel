#!/usr/bin/python2.7
import sys, getopt
import multiprocessing as mp
from ptg_utils import get_dbconnection, mean, quantile

def retrieveIG2OGscores(args):
	qogogsc, ltfamog, i, dbname, dbengine, withinfam, verbose = args
	dbco, dbcuf, dbtype, valtoken = get_dbconnection(dbname, dbengine)
	tvt = (valtoken,)*8
	q = qogogsc%tvt
	ltretval = []
	tfamog1 = ltfamog[i]
	for tfamog2 in ltfamog[i+1:]:
		if verbose: print tfamog1+tfamog2
		if tfamog1[0]==tfamog2[0]:
			if withinfam: 
				if tfamog1[1]==tfamog2[1]:
					continue
			else:
				continue
		## fetch the coevolution scores between two fam_OG groups
		dbcuf.execute(qogogsc, (tfamog1+tfamog2)*2)
		ltscores = dbcuf.fetchall()
		llt = len(ltscores)
		if llt==0:
			continue
		## filter lineage-wise bidirectional best hits
		
		# implementation 1: for each OG separately, register all lineages' best hit links
		# and then take instersection of links to ensure link is a best hit in both collection
		#~ # get best hit set for each OG as query
		#~ for k in (0,1):
			#~ if verbose: print "query OG:", [tfamog1, tfamog2][k]
			#~ stlbh = set([])
			#~ # order the results by rlocdsid of the query OG
			#~ ltscores.sort(key=lambda x: x[k])
			#~ rlocdsid = None; lbh = 0; ltlbh = []
			#~ for tscore in ltscores:
				#~ if tscore[k]!=rlocdsid:
					#~ if verbose: print "rlocdsid, ltlbh:", rlocdsid, ltlbh
					#~ # store best value for previous lineage
					#~ stlbh |= set(ltlbh)
					#~ # initiate search for next lineage
					#~ rlocdsid = tscore[k]; lbh = 0; ltlbh = []
				#~ coev = float(tscore[2])
				#~ if coev > lbh:
					#~ lbh = coev; ltlbh = [tscore]
				#~ elif coev == lbh:
					#~ # record ties for comparison with set with other OG as query
					#~ ltlbh.append(tscore)
			#~ stlbh |= set(ltlbh)
			#~ if verbose: print "stlbh (%d):"%k, stlbh
			#~ lstlbh.append(stlbh)
		#~ # merge each best-hit set
		#~ listlbh = list(lstlbh[0] & lstlbh[1])
		#~ llbh = [float(t[2]) for t in listlbh]
		#~ # this is biased by excess of many-to-many links with identical low scores
		
		# implementation 2: keep each lineage's best hit
		# unbiased but do not ensure bidirectionality; however bidirectional hits are given double weight
		dlbh = {}
		#~ dtlbh = {}
		for tscore in ltscores:
			coev = float(tscore[2])
			for k in (0,1):
				rlocdsid = tscore[k]
				prevscore = dlbh.get(rlocdsid, 0.0)
				if coev > prevscore:
					dlbh[rlocdsid] = coev
					#~ dtlbh[rlocdsid] = tscore
		# could exploit the information of who is best linked with whom
		#~ listlbh = dtlbh.values()
		#~ if verbose: print "listlbh:", listlbh
		#~ llbh = [float(t[2]) for t in listlbh]
		# but for the moment only collect scores
		llbh = dlbh.values()
		if verbose: print "llbh:", llbh
		# collect statistics
		ltretval.append(tfamog1+tfamog2+(llt, len(llbh), mean(llbh))+quantile(llbh, [0, 0.25, 0.5, 0.75, 1]))
	dbco.close()
	return ltretval

def main(orthocolid, reccolid, nfout, dbname, dbengine='postgres', withinfam=False, restrictfamogq='', nffamogqlist=None, nbthreads=1, verbose=False):
	
	# open DB connection
	dbcon, dbcur, dbtype, valtoken = get_dbconnection(dbname, dbengine)
	# create view
	dbcur.execute("drop view if exists rlocsd2og;")
	dbcur.execute("""
	create view rlocsd2og as select rlocds_id, gene_family_id, og_id
	from orthologous_groups
	inner join replacement_label_or_cds_code2gene_families using (gene_family_id, replacement_label_or_cds_code)
	where ortholog_col_id=1;
	""")
	# commit for visibility by other connections
	dbcon.commit()
	
	# fetch list of subject (fam, og) tuples
	qfamogtup = "select distinct gene_family_id, og_id from orthologous_groups "
	qfamogtup += restrictfamogq
	qfamogtup += "order by gene_family_id, og_id ;"
	dbcur.execute(qfamogtup)
	ltfamog = dbcur.fetchall()
	if nffamogqlist:
		# subset of query (fam, og) tuples
		with open(nffamogqlist, 'r') as ffamogqlist:
			ltfamogq = [tuple(line.rstrip('\n').split('\t')) for line in ffamogqlist]
			ltfamogqi = [ltfamog.index(t) for t in ltfamogq]
	else:
		# query set same as subject set
		ltfamogqi = range(len(ltfamog)-1)
	# close connection
	dbcon.close()

	if reccolid==0:
		# skip check multiplicity of reconciliation collection
		wrc = ";"
	else:
		# check multiplicity of reconciliation collection
		dbcur.execute("select distinct reconciliation_id from phylogeny.coevolution_scores;")
		lreccol = dbcur.fetchall()
		if len(lreccol)==1:
			# will not avoid further overhead of query constraint on reconciliation_id
			wrc = ";"
		else:
			# will filter queries based on reconciliation collection
			wrc = "and reconciliation_id=%d ;"%reccolid

	qogogsc = """
	select least(rlocds_id_1, rlocds_id_2), greatest(rlocds_id_1, rlocds_id_2), coev_score 
	from coevolution_scores 
	inner join rlocsd2og as og1 on rlocds_id_1=og1.rlocds_id 
	inner join rlocsd2og as og2 on rlocds_id_2=og2.rlocds_id 
	where (og1.gene_family_id=%s and og1.og_id=%s and og2.gene_family_id=%s and og2.og_id=%s) 
	or (og2.gene_family_id=%s and og2.og_id=%s and og1.gene_family_id=%s and og1.og_id=%s) 
	"""+wrc

	# run the queries in parallel
	pool = mp.Pool(processes=nbthreads)
	iterargs = ((qogogsc, ltfamogq, ltfamogs, i, dbname, dbengine, withinfam, verbose) for i in ltfamogqi)
	iterlt = pool.imap_unordered(retrieveIG2OGscores, iterargs, chunksize=1)
	# an iterator is returned by imap_unordered()
	# one needs to actually iterate over it to have the pool of parrallel workers to compute
	fout = open(nfout, 'w')
	for lt in iterlt:
		for t in lt:
			fout.write('\t'.join([str(f) for f in t])+'\n')

	fout.close()
	
def usage():
	s = "Usage: [HELP MESSAGE INCOMPLETE]\n"
	s += "python %s {--postgresql_db dbname | --sqlite_db dbfile} --out tablefiledest [OTHER OPTIONS]\n"%sys.argv[0]
	return s

if __name__=='__main__':

	opts, args = getopt.getopt(sys.argv[1:], 'hv', ['ortho_col_id=', 'reconciliation_id=', 'out=', \
	                                                'postgresql_db=', 'sqlite_db=', 'whitinfam', \
	                                                'restrict_famog_query=', 'input_famog_query_list=', \
	                                                'threads=', 'help', 'verbose'])
	dopt = dict(opts)
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)

	nfout = dopt['--out']
	dbname = dopt.get('--postgresql_db', dopt.get('--sqlite_db'))
	if '--postgresql_db' in dopt:
		dbengine = 'postgres'
	elif '--sqlite_db' in dopt:
		dbengine = 'sqlite'
		raise NotImplementedError, "condensation of coevolution network from gene-lineage to orthologous groups not implement with SQLite database - and unlikely to be so given the huge size of tables to be loaded in memory"
	else:
		raise ValueError, "must provide database name (postgreSQL) / file location (SQLite) through '--sqlite_db' or '--postgresql_db' options"

	orthocolid = int(dopt.get('--ortho_col_id', 1))
	reccolid = int(dopt.get('--reconciliation_id', 0))
	withinfam = bool(int(dopt.get('--whitinfam', 0)))
	restrictfamogq = dopt.get('--restrict_famog_query', '')
	nffamogqlist = dopt.get('--input_famog_query_list')
	nbthreads = int(dopt.get('--threads', 1))
	verbose = (('-v' in dopt) or ('--verbose' in dopt))
	if verbose: print "dopt:", dopt
	
	main(orthocolid, reccolid, nfout, dbname, dbengine, withinfam, restrictfamogq, nffamogqlist, nbthreads, verbose)
