#!/usr/bin/python
import glob, os, sys
import psycopg2

sqldbname = sys.argv[1]
dirortho = sys.argv[2]
orthomethod = sys.argv[3]
clustmethod = sys.argv[4]
ortcolid = int(sys.argv[5])

dbcon = psycopg2.connect("dbname=%s"%sqldbname)
dbcur = dbcon.cursor()
dbcur.execute("DELETE FROM phylogeny.orthologous_groups WHERE ortholog_col_id=%s;", (ortcolid,))
filesuffix = "_%s.orthologs.%s"%(orthomethod, clustmethod)
lnfortho = glob.glob(os.path.join(dirortho, orthomethod, '*'+filesuffix))
for nfortho in lnfortho:
  fam = os.path.basename(nfortho).replace(filesuffix, '')
  with open(nfortho, 'r') as fortho:
    lcdsog = [tuple(line.replace(' ', '').rstrip('\n').split('\t')) for line in fortho]
  dbcur.executemany("INSERT INTO phylogeny.orthologous_groups (replacement_label_or_cds_code, gene_family_id, og_id, ortholog_col_id) VALUES (%s,%s,%s,%s);", [(cds, fam, ogid, ortcolid) for cds, ogid in lcdsog])

dbcur.execute("CREATE INDEX IF NOT EXISTS og_cds_idx ON phylogeny.orthologous_groups (replacement_label_or_cds_code)")
dbcur.execute("CREATE INDEX IF NOT EXISTS og_fam_idx ON phylogeny.orthologous_groups (gene_family_id)")
dbcur.execute("CREATE INDEX IF NOT EXISTS og_fam_ogid_idx ON phylogeny.orthologous_groups (gene_family_id, og_id)")
dbcur.execute("CREATE UNIQUE INDEX IF NOT EXISTS og_cds_ogcol_idx ON phylogeny.orthologous_groups (replacement_label_or_cds_code, ortholog_col_id)")

dbcon.commit()
dbcon.close()
