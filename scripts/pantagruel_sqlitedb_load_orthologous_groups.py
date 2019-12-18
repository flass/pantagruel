#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import glob, os, sys
import sqlite3

nfsqldb = sys.argv[1]
dirortho = sys.argv[2]
orthomethod = sys.argv[3]
clustmethod = sys.argv[4]
ortcolid = int(sys.argv[5])
if len(sys.argv)>6:
	keepPrevRecord = bool(int(sys.argv[6]))
else:
	keepPrevRecord = False

dbcon = sqlite3.connect(nfsqldb)
dbcur = dbcon.cursor()
if not keepPrevRecord: dbcur.execute("DELETE FROM orthologous_groups WHERE ortholog_col_id=?;", (ortcolid,))
filesuffix = "_%s.orthologs.%s"%(orthomethod, clustmethod)
lnfortho = glob.glob(os.path.join(dirortho, orthomethod, '*'+filesuffix))
for nfortho in lnfortho:
    fam = os.path.basename(nfortho).replace(filesuffix, '')
    with open(nfortho, 'r') as fortho:
        lcdsog = [tuple(line.replace(' ', '').rstrip('\n').split('\t')) for line in fortho]
    dbcur.executemany("INSERT INTO orthologous_groups (replacement_label_or_cds_code, gene_family_id, og_id, ortholog_col_id) VALUES (?,?,?,?);", [(cds, fam, ogid, ortcolid) for cds, ogid in lcdsog])
    #~ dbcur.executemany("INSERT INTO orthologous_groups (cds_code, gene_family_id, og_id, ortholog_col_id) VALUES (?,?,?,?);", [(cds, fam, ogid, ortcolid) for cds, ogid in lcdsog])

dbcur.execute("CREATE INDEX IF NOT EXISTS og_cds_idx ON orthologous_groups (replacement_label_or_cds_code);")
#~ dbcur.execute("CREATE INDEX IF NOT EXISTS og_cds_idx ON orthologous_groups (cds_code);")
dbcur.execute("CREATE INDEX IF NOT EXISTS og_fam_idx ON orthologous_groups (gene_family_id);")
dbcur.execute("CREATE INDEX IF NOT EXISTS og_fam_ogid_idx ON orthologous_groups (gene_family_id, og_id);")
dbcur.execute("CREATE UNIQUE INDEX IF NOT EXISTS og_cds_ogcol_idx ON orthologous_groups (replacement_label_or_cds_code, ortholog_col_id);")
#~ dbcur.execute("CREATE UNIQUE INDEX IF NOT EXISTS og_cds_ogcol_idx ON orthologous_groups (cds_code, ortholog_col_id);")

dbcur.execute("DROP TABLE IF EXISTS gene_fam_og_sizes;")
dbcur.execute("""CREATE TABLE gene_fam_og_sizes AS 
                   SELECT gene_family_id, og_id, count(cds_code) AS size, count(cds_code) AS genome_present, ortholog_col_id 
                    FROM orthologous_groups 
                    INNER JOIN genome.gene_tree_label2cds_code USING (replacement_label_or_cds_code)
                   WHERE ortholog_col_id=?
                   GROUP BY gene_family_id, og_id, ortholog_col_id
                  UNION
                    SELECT gene_family_id, CAST(NULL AS INT) AS og_id, size, genome_present, ? AS ortholog_col_id 
                     FROM gene_family_sizes
                    WHERE gene_family_id NOT IN (SELECT DISTINCT gene_family_id FROM orthologous_groups)
                  ;""", (ortcolid, ortcolid))

#~ dbcur.execute("""CREATE TABLE gene_fam_og_sizes AS 
                   #~ SELECT gene_family_id, og_id, count(cds_code) AS size, count(cds_code) AS genome_present, ortholog_col_id 
                    #~ FROM orthologous_groups 
                   #~ WHERE ortholog_col_id=%d
                   #~ GROUP BY gene_family_id, og_id, ortholog_col_id
                  #~ UNION
                    #~ SELECT gene_family_id, CAST(NULL AS INT) AS og_id, size, genome_present, %d AS ortholog_col_id 
                     #~ FROM gene_family_sizes
                    #~ WHERE gene_family_id NOT IN (SELECT DISTINCT gene_family_id FROM orthologous_groups)
                  #~ ;"""%(ortcolid, ortcolid))

dbcur.execute("CREATE INDEX IF NOT EXISTS og_size_size_idx ON gene_fam_og_sizes (size);")
dbcur.execute("CREATE INDEX IF NOT EXISTS og_size_present_idx ON gene_fam_og_sizes (genome_present);")
dbcur.execute("CREATE INDEX IF NOT EXISTS og_size_famog_idx ON gene_fam_og_sizes (gene_family_id, og_id);")

dbcon.commit()
dbcon.close()
