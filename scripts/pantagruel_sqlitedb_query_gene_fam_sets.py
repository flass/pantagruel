#/usr/bin/python
import sqlite3
import sys, os, getopt, re

def main(dbconpar, baseq, minsizes=[], maxsizes=[], dirout='', outprefix='', nfexcludefams=None):
	
	assert len(minsizes)==len(maxsizes)
	dbcon = sqlite3.connect(dbconpar)
	dbcur = dbcon.cursor()
	
	if nfexcludefams:
		dbcur.execute( 'create temporary table excludefams (gene_family_id varchar(25));' )
		with open(nfexcludefams, 'r') as fexcludefams:
			cursor.executemany("INSERT INTO excludefams VALUES (?);", ((line.rstrip('\n'),) for line in fexcludefams))
		joinexcl = " LEFT JOIN excludefams USING (gene_family_id)"
		whereexcl = " AND excludefams.gene_family_id IS NULL"
		basqparts = re.split(' where ', baseq.rstrip(';'), maxsplit=1, flags=re.IGNORECASE)
		headq = basqparts[0]+joinexcl+' WHERE '+basqparts[1]+whereexcl
	else:
		headq = baseq.rstrip(';')
	
	for i in range(len(minsizes)):
		try:
			mins = int(minsizes[i])
			gtmin = " and size>=%d"%mins
			smins = "_minsize%d"%mins
		except ValueError:
			gtmin = ""
			smins = ""
		try:
			maxs = int(maxsizes[i])
			ltmax = " and size<=%d"%maxs
			smaxs = "_maxsize%d"%maxs
		except ValueError:
			ltmax = ""
			smaxs = ""
		
		q = "%s%s%s order by size desc;"%(baseq, gtmin, ltmax)
		dbcur.execute( q )
		lfamsizes = dbcur.fetchall()
		with open(os.path.join(dirout, "%s%s%s"%(outprefix, smins, smaxs)), 'w') as fout:
			for row in lfamsizes:
				fout.write('\t'.join([str(x).rstrip(' ') for x in row])+'\n')

def usage():
	s = 'Usage: python %s --db=sqlitedb_file_path [ --famsets.min.sizes="[x][,y,...,]z" --famsets.max.sizes=X[,Y,...,][Z] ] [OPTIONS]'%(sys.argv[0])
	return s
	
if __name__=='__main__':	
	opts, args = getopt.getopt(sys.argv[1:], 'h', ['famsets.min.sizes=', 'famsets.max.sizes=', \
	                                               'db=', 'exclude.fams=', \
	                                               'base.query=', 'cds.orfan.fam.name=', \
	                                               'dirout=', 'outprefix=', \
	                                               'help'])
	
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	db = dopt['--db']
	minsizes = dopt.get('--famsets.min.sizes', '').split(',')
	maxsizes = dopt.get('--famsets.max.sizes', '').split(',')
	cdsorfanclust = dopt.get('--cds.orfan.fam.name', '')
	basequery = dopt.get('--base.query', "select gene_family_id, size from gene_family_sizes where gene_family_sizes.gene_family_id is not null and gene_family_sizes.gene_family_id!='%s'"%(cdsorfanclust))
	dirout = dopt.get('--dirout', os.getcwd())
	outprefix = dopt.get('--outprefix', 'cdsfam')
	nfexcludefams = dopt.get('--exclude.fams')
	
	main(db, baseq=basequery, minsizes=minsizes, maxsizes=maxsizes, dirout=dirout, outprefix=outprefix, nfexcludefams=nfexcludefams)
	
