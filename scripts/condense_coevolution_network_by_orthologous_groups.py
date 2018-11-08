#!/usr/bin/python2.7
import psycopg2
import sys
import multiprocessing as mp
from ptg_utils import get_dbconnection, mean, median

sqldbname = sys.argv[1]
orthocolid = int(sys.argv[2])
reccolid = int(sys.argv[3])
nfout = sys.argv[4]
if len(sys.argv)>5:
	withinfam = bool(int(sys.argv[5]))
else:
	withinfam = False

dbcon = psycopg2.connect("dbname=%s"%sqldbname)
dbcur = dbcon.cursor()

dbcur.execute("drop view if exists rlocsd2og;")
dbcur.execute("""
create view  as select rlocds_id, gene_family_id, og_id
from orthologous_groups
inner join replacement_label_or_cds_code2gene_families using (gene_family_id, replacement_label_or_cds_code)
where ortholog_col_id=1;
""")

dbcur.execute("select distinct gene_family_id, og_id from orthologous_groups order by gene_family_id, og_id")
ltfamog = dbcur.fetchall()

# check multiplicity of reconciliation collection
dbcur.execute("select distinct reconciliation_id from phylogeny.coevolution_scores;")
lreccol = dbcur.fetchall()
if len(lreccol)==1:
	# will avoid further overhead of query constraint on reconciliation_id
	wrc = ";"
else:
	wrc = "and reconciliation_id=%d ;"%reccolid

qogogsc = """select rlocds_id_1, rlocds_id_2, coev_score 
from coevolution_scores 
inner join rlocsd2ogcol1 as og1 on rlocds_id_1=og1.rlocds_id 
inner join rlocsd2ogcol1 as og2 on rlocds_id_2=og2.rlocds_id 
where  og1.gene_family_id=%s and og1.rlocds_id=%s and og2.gene_family_id=%s and og2.rlocds_id=%s 
"""+wrc

def retrieveIG2OGscores(qogogsc, ltfamog, i):
	dbcon, dbcuf, dbtype, valtoken = get_dbconnection(dbname, dbengine)
	q = qogogsc%(valtoken, valtoken, valtoken, valtoken)
	lretval = []
	tfamog1 = ltfamog[i]
	for tfamog2 in ltfamog[i+1:]:
		if tfamog1[0]==tfamog2[0]:
			if withinfam: 
				if tfamog1[1]==tfamog2[1]:
					continue
			else:
				continue
		## fetch the coevolution scores between two fam_OG groups
		dbcuf.execute(qogogsc, tfamog1+tfamog2)
		ltscores = dbcuf.fetchall()
		if len(ltscores)==0:
			continue
		## filter lineage-wise bidirectional best hits
		lstlbh = []
		# get best hit set for each OG as query
		for k in (0,1):
			stlbh = set([])
			# order the results by rlocdsid of the query OG
			ltscores.sort(key=lambda x: x[k])
			rlocdsid = tscore[k]; lbh = 0; ltlbh = []
			for tscore in ltscores:
				if tscore[k]!=rlocdsid:
					# store best value for previous lineage
					stlbh |= set(ltlbh)
					# initiate search for next lineage
					rlocdsid = tscore[k]; lbh = 0; ltlbh = []
				coev = float(tscore[2])
				if coev > lbh:
					lbh = coev; ltlbh = [tscore]
				elif coev == lbh:
					# record ties for comparison with set with other OG as query
					ltlbh.append(tscore)
		# merge each best-hit set
		stlbh = lstlbh[0] & lstlbh[1]
		# collect statistics
		llbh = [float(t[2]) for t in stlbh]
		# could exploit the information of who is best linked with whom, but not for the moment
		minlbh = min(llbh)
		maxlbh = max(llbh)
		meanlbh = mean(llbh)
		medlbh = median(llbh)
		ltretval.append(tfamog1+tfamog2+(len(llbh), minlbh, maxlbh, meanlbh, medlbh)))
	dbcon.close()
	return ltretval

# run the queries in parallel
pool = mp.Pool(processes=nbthreads)
iterargs = ((qogogsc, ltfamog, i) for i in len(ltfamog)-1)
iterlt = pool.imap_unordered(retrieveIG2OGscores, iterargs, chunksize=1)
# an iterator is returned by imap_unordered()
# one needs to actually iterate over it to have the pool of parrallel workers to compute
fout = open(nfout, 'w')
for ltretval in iterlt:
	for t in lt:
		fout.write('\t'.join([str(f) for f in t])+'\n')
