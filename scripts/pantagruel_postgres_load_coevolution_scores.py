#!/usr/bin/python
import glob, os, sys
import psycopg2

fileprefix = 'matching_events.tab'

sqldbname = sys.argv[1]
dirscores = sys.argv[2]
reccolid = int(sys.argv[3])

dbcon = psycopg2.connect("dbname=%s"%sqldbname)
dbcur = dbcon.cursor()
dbcur.execute("DELETE FROM phylogeny.coevolution_scores WHERE reconciliation_id=%s;", (reccolid,))
nfscorepat = os.path.join(dirscores, fileprefix+'*')
lnfscore = glob.glob(nfscorepat)
for nfscore in lnfscore:
  print nfscore
  with open(nfscore, 'r') as fscore:
	dbcur.copy_from(file=fscore, table='phylogeny.coevolution_scores', sep='\t', size=65536, columns=('rlocds_id_1', 'rlocds_id_2', 'coev_score'))

dbcon.commit()
print "COMMITED COPY phylogeny.coevolution_scores FROM %s"%nfscorepat

dbcur.execute("CREATE UNIQUE INDEX IF NOT EXISTS rlocds_id_pair_idx ON phylogeny.coevolution_scores (rlocds_id_1, rlocds_id_2)")
dbcon.commit()
print "COMMITED CREATE UNIQUE INDEX rlocds_id_pair_idx"

dbcur.execute("UPDATE phylogeny.coevolution_scores SET reconciliation_id=%s WHERE reconciliation_id IS NULL;", (reccolid,))
dbcur.execute("CREATE INDEX IF NOT EXISTS reccolid_idx ON phylogeny.coevolution_scores (reconciliation_id)")
dbcon.commit()
print "COMMITED UPDATE SET reconciliation_id=%d AND CREATE INDEX reccolid_idx"%reccolid

dbcon.close()
