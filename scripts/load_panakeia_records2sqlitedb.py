#!/usr/bin/env python2.7
import sys, os, getopt
import sqlite3
import time

def main(nfsqldb, nfprotclust, nfprotsubclust=None, nfsynpat=None, subclusogcolid=None, verbose=False):
	# prepare database
	dbcon = sqlite3.connect(nfsqldb, isolation_level="DEFERRED")
	dbcur = dbcon.cursor()
	if subclusogcolid:
		# first check the database is clear for that collection id
		dbcur.execute("SELECT ortholog_col_id FROM ortholog_collections WHERE ortholog_col_id=?;", subclusogcolid)
		ogcolidexists = dbcur.fetchall()
		if ogcolidexists:
			raise IndexError, "orthologous collection with id %d is already used in the database '%s':\n%s\n- - -\nPlease use another integer as input OG collection id."%(subclusogcolid, nfsqldb, '\n'.join([repr(row) for row in ogcolidexists])')
	
	# parse input data
	dprotclust = {}
	dprotsubclust = {}
	dsynpat = {}
	with open(nfprotclust, 'r') as fprotclust:
		for line in fprotclust:
			lsp = line.rstrip('\n').split(': ')
			dprotclust[lsp[0]] = lsp[1].split('\t')
	lprotclust = sorted(dprotclust.keys())
	if nfprotsubclust:
		with open(nfprotsubclust, 'r') as fprotsubclust:
			for line in fprotsubclust:
				lsp = line.rstrip('\n').split(': ')
				dprotsubclust[lsp[0]] = lsp[1].split('\t')
		lprotsubclust = sorted(dprotsubclust.keys())
	if nfsynpat:
		with open(nfsynpat, 'r') as fsynpat:
			for line in fsynpat:
				lsp = line.rstrip('\n').split(': ')
				patit, pattype = lsp[0].split(', ')	# TO ADAPT TO ACTUAL INPUT FORMAT
				dsynpat[(patit, pattype)] = lsp[1].split('\t')
		lsynpat = sorted(dsynpat.keys())
	
	# (re)set tables
	dbcur.execute("DROP TABLE IF EXISTS panakeia_gene_clusters;")
	dbcur.execute("DROP TABLE IF EXISTS panakeia_gene_subclusters;")
	dbcur.execute("DROP TABLE IF EXISTS panakeia_gene_patterns;")
	dbcur.execute("""CREATE TABLE panakeia_gene_clusters (
	                   cds_code varchar(20),
					   cluster_id varchar(20)
	                 );""") 
	dbcur.execute("""CREATE TABLE panakeia_gene_subclusters (
	                   cds_code varchar(20),
					   subcluster_id varchar(20)
	                 );""") 
	dbcur.execute("""CREATE TABLE panakeia_gene_patterns (
	                   cds_code varchar(20),
					   pattern_id varchar(20)
					   pattern_type varchar(20)
	                 );""")
	# load table contents
	dbcur.executemany("""INSERT INTO panakeia_gene_clusters VALUES (?,?)""", ((cdscode, clusid) for clusid in lprotclust for cdscode in dprotclust[clusid]))
	if nfprotsubclust:
		dbcur.executemany("""INSERT INTO panakeia_gene_subclusters VALUES (?,?)""", ((cdscode, sclusid) for sclusid in lprotsubclust for cdscode in dprotsubclust[sclusid]))
	if nfsynpat:
		dbcur.executemany("""INSERT INTO panakeia_gene_patterns VALUES (?,?,?)""", ((cdscode, patid, pattype) for patid,pattype in lsynpat for cdscode in dsynpat[(patid,pattype)]))
	
	if subclusogcolid:
		# first set the new OG collection
		dbcur.execute("INSERT INTO ortholog_collections VALUES (?,?,?,?,?,?,?,datetime('now'));", (subclusogcolid, 'panakeia_subclusters', 0, 'Panakeia', 'v0.1', 'main'))
		
		# query the mapping of CDS to Pantagruel gene families
		qfam = "SELECT gene_family_id FROM coding_sequences;"
		dbcur.execute(qfam)
		lfam = dbcur.fetchall()
		for fam in lfam:
			# retrieve CDS and associated Panakeia subclusters
			if nfsynpat:
				qcdssubc = "SELECT cds_code, pattern_id FROM coding_sequences LEFT JOIN panakeia_gene_patterns USING (cds_code) WHERE gene_family_id=? AND pattern_type='indel';"
			elif nfprotsubclust:
				qcdssubc = "SELECT cds_code, subcluster_id FROM coding_sequences LEFT JOIN panakeia_gene_subclusters USING (cds_code) WHERE gene_family_id=?;"
			else:
				qcdssubc = "SELECT cds_code, cluster_id FROM coding_sequences LEFT JOIN panakeia_gene_clusters USING (cds_code) WHERE gene_family_id=?;"
			dbcur.execute(qcdssubc, fam)
			lcdssubc = dbcur.fetchall()
			usubc = set([cdssubc[1] for cdssubc in lcdssubc])
			dsubcid = dict((subc, i) for i, subc in enumerate(usubc))
			dbcur.executemany("INSERT INTO orthologous_groups VALUES (?,?,?,?)", ((cds, fam, dsubcid[subc], subclusogcolid) for cds,subc in lcdssubc))
	
	dbcon.commit()

def usage():
	s =  'Basic usage:\n'
	s += ' python %s -d /path/to/sqlite3.db.file -p /path/to/protein.cluster.file -s /path/to/protein.subcluster.file -t /path/to/synteny.pattern.file [options]\n'%sys.argv[0]
	s += 'Options:\n'
	s += '  -i|--otholog_col_id int\tspecify id of orthologous group collection to which Panakeia subclusters will be assigned.\n'
	s += '  -v|--verbose\t\t\tverbose mode.\n'
	s += '  -h|--help\t\t\tthis message.\n'
	return s


if __name__=='__main__':
	
	opts, args = getopt.getopt(sys.argv[1:], 'd:p:s:t:i:vh', ['verbose', 'help', 'otholog_col_id='])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)

	dbname = dopt['-d']
	nfprotclust = dopt['-p']
	nfprotsubclust = dopt['-s']
	nfsynpat = dopt['-t']
	subclusogcolid = dopt.get('-i', dopt.get('--otholog_col_id'))
	verbose = (('--verbose' in dopt) or ('-v' in dopt))
	
	for nf in [dbname, nfprotclust, nfprotsubclust, nfsynpat]:
		if (nf is not None) and not (os.path.exists(nf)):
			raise ValueError, "specified input file '%s' cannot be found"%nf
	
	main(dbname, nfprotclust, nfprotsubclust, nfsynpat, subclusogcolid=subclusogcolid, verbose=verbose)
