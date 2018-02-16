#!/usr/bin/python
import re
import os, sys
import sqlite3

def replaceValuesAsNull(table, cursor, nullval='', ommitcols=[]):
	cursor.execute("PRAGMA table_info(%s);"%table)
	colinfos = cursor.fetchall()
	for colinfo in colinfos:
		colname = colinfo['name']
		cursor.execute("UPDATE %s SET %s=NULL WHERE %s='%s';"%(table, colname, colname, nullval)

def loadAndCurateTable(table, nfin, cursor, insertcolumns=(), sep='\t', **kw):
	if insertcolumns:
		if type(insertcolumns)==str: coldef = insertcolumns
		else: coldef = "(%s)"%(', '.join(insertcolumns))
	else:
		coldef = ''
	with open(nfin, 'r') as ftabin:
		cursor.executemany("INSERT INTO %s %s VALUES (?, ?);"%(table, coldef), tuple(line.split(sep) for line in ftabin))
	replaceValuesAsNull(table, cursor, **kw)
	
def createAndLoadTable(table, tabledef, nfin, cursor, temp=False, enddrop=False, **kw):
	tmp = 'TEMP' if temp else ''
	cursor.execute("CREATE %s TABLE %s %s ON COMMIT DROP %;"%(tmp, table, tabledef, tmpdrop))
	loadAndCurateTable(table, nfin, cursor, **kw)
	if enddrop: cursor.execute("DROP TABLE %s;"%table)

dbname = sys.argv[1]
protfamseqs = os.environ['protfamseqs']
protorfanclust = os.environ['protorfanclust']
cdsorfanclust = os.environ['cdsorfanclust']

conn = sqlite3.connect(database=dbname)
cur = conn.cursor()


# populate assemblies table
loadAndCurateTable('assemblies', 'genome_assemblies.tab', cur, ommitcols=['assembly_id', 'assembly_name'])
# populate replicons table
loadAndCurateTable('replicons', 'genome_replicons.tab', cur, insertcolumns="(assembly_id, genomic_accession, replicon_name, replicon_type, replicon_size)")
conn.commit()

# populate protein table
protprodtabledef = """(genbank_nr_protein_id VARCHAR(20), product TEXT)"""
protfamtabledef = """(genbank_nr_protein_id VARCHAR(20), protein_family_id CHAR(13))"""
createAndLoadTable('protein_products', protprodtabledef, 'genome_protein_products.tab', cur)
createAndLoadTable('protein_fams', protfamtabledef, protfamseqs, cur)
cur.executescript("""
INSERT INTO proteins (genbank_nr_protein_id, product, protein_family_id) (select genbank_nr_protein_id, min(product), protein_family_id
 FROM protein_products
 INNER JOIN protein_fams USING (genbank_nr_protein_id)
 GROUP BY genbank_nr_protein_id,protein_family_id
);
DROP TABLE protein_products;
DROP TABLE protein_fams;
""")
conn.commit()

# populate coding_sequences table
cdstabledef = """(
  nr_protein_id CHAR(15),
  genomic_accession CHAR(14) NOT NULL,
  locus_tag VARCHAR(200),
  cds_begin INTEGER NOT NULL,
  cds_end INTEGER NOT NULL,
  strand CHAR(1) NOT NULL,
  genbank_cds_id VARCHAR(50) NOT NULL
)"""
cdsfamtabledef = """(
  genbank_cds_id VARCHAR(50) NOT NULL,
  gene_family_id CHAR(13)
)"""
createAndLoadTable('codingsequences', cdstabledef, 'genome_coding_sequences.tab', cur)
createAndLoadTable('cdsfam', cdsfamtabledef, 'genome_gene_families.tab', cur)
cur.executescript("""
CREATE INDEX gbcdsid_1 ON codingsequences (genbank_cds_id);
CREATE INDEX gbcdsid_2 ON cdsfam (genbank_cds_id);
INSERT INTO coding_sequences (genbank_cds_id, genomic_accession, locus_tag, cds_begin, cds_end, strand, genbank_nr_protein_id, gene_family_id) (
 SELECT genbank_cds_id, genomic_accession, locus_tag, cds_begin, cds_end, strand, nr_protein_id, gene_family_id
  FROM codingsequences
  LEFT JOIN cdsfam USING (genbank_cds_id)
);
"""
conn.commit()

# populate the protein families
cur.execute("INSERT INTO nr_protein_families (SELECT DISTINCT protein_family_id, false FROM proteins);")
# add the bit mark for the singleton nr protein family
cur.execute("UPDATE nr_protein_families SET is_singleton=true WHERE protein_family_id=?;", protorfanclust)
# allocated the '$cdsorfanclust' value to gene_family_id field for those CDSs with a parent protein but no gene family affiliation
cur.execute("UPDATE coding_sequences SET gene_family_id=? " \
           +"WHERE gene_family_id IS NULL AND genbank_nr_protein_id IS NOT NULL;", cdsorfanclust)
# populate the gene families deriving from the protein families
cur.execute("INSERT INTO gene_families " \
           +"(SELECT replace(protein_family_id, 'P', 'C') AS gene_family_id, false, protein_family_id " \
           +"FROM nr_protein_families);")
# add the gene families which have no natural parent nr protein family, i.e. those derived from singleton nr proteins, but with multiple CDS members
cur.execute("INSERT INTO gene_families (SELECT distinct gene_family_id, false, ? " \
           +"FROM coding_sequences LEFT JOIN gene_families USING (gene_family_id) " \
           +"WHERE protein_family_id IS NULL AND gene_family_id IS NOT NULL);", protorfanclust)
# add the bit mark for the ORFan gene family
cur.execute("UPDATE gene_families SET is_orfan=true WHERE gene_family_id=?;", cdsorfanclust
# finally set the primary key
cur.execute("ALTER TABLE gene_families ADD PRIMARY KEY (gene_family_id);"
