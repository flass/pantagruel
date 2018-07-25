#!/usr/bin/python
import glob, os, sys
import sqlite3

nfsqldb = sys.argv[1]
dirortho = sys.argv[2]
orthomethod = sys.argv[3]
clustmethod = sys.argv[4]
ortcolid = int(sys.argv[5])

dbcon = sqlite3.connect(nfsqldb)
dbcur = dbcon.cursor()
dbcur.execute("DELETE FROM orthologous_groups WHERE ortholog_col_id=?;", (ortcolid,))
filesuffix = "_%s.orthologs.%s"%(orthomethod, clustmethod)
lnfortho = glob.glob(os.path.join(dirortho, orthomethod, '*'+filesuffix))
for nfortho in lnfortho:
  fam = os.path.basename(nfortho).replace(filesuffix, '')
  with open(nfortho, 'r') as fortho:
    lcdsog = [tuple(line.replace(' ', '').rstrip('\n').split('\t')) for line in fortho]
  dbcur.executemany("INSERT INTO orthologous_groups (cds_code, gene_family_id, og_id, ortholog_col_id) VALUES (?,?,?,?);", [(cds, fam, ogid, ortcolid) for cds, ogid in lcdsog])

dbcur.execute("CREATE INDEX og_cds_idx ON orthologous_groups (cds_code)")
dbcur.execute("CREATE INDEX og_fam_idx ON orthologous_groups (gene_family_id)")
dbcur.execute("CREATE INDEX og_fam_ogid_idx ON orthologous_groups (gene_family_id, og_id)")
dbcur.execute("CREATE UNIQUE INDEX og_cds_ogcol_idx ON orthologous_groups (cds_code, ortholog_col_id)")

dbcon.commit()
dbcon.close()
