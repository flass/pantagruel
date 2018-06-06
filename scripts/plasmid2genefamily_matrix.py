#!/usr/bin/python
import psycopg2
import os, sys, getopt

def main(nfout, minoccur, labexcl, dbname, dbuser, dbhost):
	
	conn = psycopg2.connect(dbname=dbname, user=dbuser, host=dbhost)
	cur = conn.cursor()

	fout = open(nfout, 'w')
	if exclpat:
		foutexcl = open("%s_excl-%s"%(nfout, labexcl), 'w')
	
	cur.execute("select genomic_accession from genome.replicons where replicon_type='plasmid';")
	ltplasmidacc = cur.fetchall()
	lplasmidacc = [tplasmidacc[0].rstrip(' ') for tplasmidacc in ltplasmidacc]
	lplasmidacc.sort()

	dplas2fam = {}
	for pacc in lplasmidacc:
		print pacc
		cur.execute("select gene_family_id, count(*) from genome.coding_sequences where genomic_accession=%s and gene_family_id is not null group by gene_family_id;", (pacc,))
		lplasfams = cur.fetchall()
		for tplasfam in lplasfams:
			fam, famcount = tplasfam ; fam = fam.rstrip(' ')
			dplas2fam.setdefault(fam, {}).setdefault(pacc, int(famcount))

	lfam = dplas2fam.keys()
	lfam.sort()
	if exclpat:
		cur.execute("select distinct(gene_family_id) from genome.proteins inner join genome.gene_families using (protein_family_id) where product ~* '%s' ;"%exclpat)
		ltfamexcl = cur.fetchall()
		lfamexcl = [tfamexcl[0].rstrip(' ') for tfamexcl in ltfamexcl]
	else:
		lfamexcl = []

	#header
	fout.write('\t'.join(['']+lplasmidacc)+'\n')
	if exclpat: foutexcl.write('\t'.join(['']+lplasmidacc)+'\n')
	#rows = families
	for fam in lfam:
		if len(dplas2fam[fam])>=minoccur:
			sout = '\t'.join([fam]+[str(dplas2fam[fam].get(pacc, 0)) for pacc in lplasmidacc])+'\n'
			fout.write(sout)
			if exclpat and (fam not in lfamexcl):
				foutexcl.write(sout)
		
	fout.close()
	if exclpat: foutexcl.close()		

	
def usage():
	s =  ''
	return s
	
if __name__=='__main__':	
	opts, args = getopt.getopt(sys.argv[1:], 'o:h', ['out.file=', \
	                                               'min.fam.occur=', 'exclude.pat=', 'label.exclude=', \
	                                               'sql.host.name=', 'sql.db.name=', 'sql.user.name=', \
	                                               'help'])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)

	nfout = dopt.get('-o', dopt.get('--out.file'))
	if not nfout:
		print "Error: must provide output file path through options '-o' or '--output.file'"
		print usage()
		sys.exit(1)
	
	minoccur = int(dopt.get('--min.fam.occur', 2))
	exclpat = dopt.get('--exclude.pat', '(tnp|transpos|insertion (element|sequence))')
	labexcl = dopt.get('--label.exclude', 'transposases')
	
	dbname = dopt.get('--sql.db.name', os.environ['sqldbname'])
	dbuser = dopt.get('--sql.user.name', os.environ['USER'])
	dbhost = dopt.get('--sql.host.name', 'localhost')
	
	main(nfout=nfout, minoccur=minoccur, labexcl=labexcl, dbname=dbname, dbuser=dbuser, dbhost=dbhost)
