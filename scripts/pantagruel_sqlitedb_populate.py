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
curry = conn.cursor()


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

# load UniProt taxon codes for CDS name shortening
fout = open("uniprotcode_taxid.tab", 'w')
codetaxids = []
CODepat = re.compile("([A-Z0-9]{3,5}) +[ABEVO] +([0-9]{1,7}): ")
fspeclist = open('speclist', 'r')
for line in fspeclist:
  if line[0]==' ': continue
  CODEmatch = CODepat.match(line)
  if CODEmatch:
    code, taxid = CODEmatch.groups()
    fout.write("%s\t%s\n"%(code, taxid))
    codetaxids.append((code, taxid))

fout.close()
fspeclist.close()
cur.execute("INSERT INTO uniptrotcode2taxid VALUES (?, ?);", codetaxids)

# generate unique code for each genome assembly from the Uniprot code + enough digits
cur.execute("SELECT assembly_id, uniptrotcode2taxid.code, species FROM assemblies LEFT JOIN uniptrotcode2taxid USING (taxid);")
lasscode = cur.fetchall()
try:
  cur.execute("ALTER TABLE assemblies ADD COLUMN code varchar(10);")
except sqlite3.Error, e:
  conn.rollback()
  sys.stderr.write( 'Warning:' + str(e) + "reset assemblies.code values to NULL before update\n" )
  cur.execute("UPDATE assemblies SET code=NULL;")

dcodesn = {}
for ass, code, spe in lasscode:
  if not code:
    s, p = spe.split()
    if '.' in p:
      # fairly common case of organism name being 'Genus sp.'
      c = s[:6].upper()
    else:
      # regular case  of organism name being 'Genus species'
      c = (s[:3]+p[:3]).upper()
  else:
    c = code
  dcodesn[c] = dcodesn.setdefault(c, 0) + 1
  cur.execute("update assemblies set code=? where assembly_id=?;", (c+str(dcodesn[c]), ass.strip(' ')))
  print c+str(dcodesn[c]), ass
conn.commit()

# add the cds_code column to coding_sequences table, generate unique cds codes based on the species code and genbank CDS id and create indexes
cur.execute("PRAGMA table_info(coding_sequences);")
cdscolinfos = cursor.fetchall()
for cdscoli, cdscolinfo in cdscolinfos:
	if cdscolinfo['name']=='genbank_cds_id': break
else:
	raise ValueError, "column 'genbank_cds_id' is missing in table coding_sequences"

cur.execute("""
SELECT coding_sequences.*, code
FROM coding_sequences
INNER JOIN replicons using (genomic_accession)
INNER JOIN assemblies USING (assembly_id);
""")
cdsrows = curry.fetchall()
csdnumpat = re.compile('.+_([0-9]+)$')
csdnumrows = (cdsrow[cdscoli:]+[cdsnumpat.search(cdsrow[cdscoli]).group(1)]+[:cdscoli+1] for cdsrow in cdsrows)

cur.executescript("""
create table coding_sequences2 (like coding_sequences);"
alter table coding_sequences2 add column cds_code varchar(20);
""")
cur.executemany("INSERT INTO coding_sequences2 VALUES (%s)"%(','.join(['?']*len(cdscolinfo))), csdnumrows)
cur.executescript("""
drop table coding_sequences;
alter table coding_sequences2 rename to coding_sequences;
-- recreate indexes and views
CREATE INDEX genbank_cds_id_key ON coding_sequences (genbank_cds_id);
CREATE INDEX cds_code_key ON coding_sequences (cds_code);
CREATE INDEX gene_family_id ON coding_sequences (gene_family_id);
CREATE INDEX genomic_accession_key ON coding_sequences (genomic_accession);
CREATE INDEX genomic_accession_cds_begin_key ON coding_sequences (genomic_accession, cds_begin);
CREATE INDEX genbank_nr_protein_id_key ON coding_sequences (genbank_nr_protein_id);
CREATE VIEW gene_family_sizes AS SELECT gene_family_id, count(*) as size from coding_sequences group by gene_family_id;
ALTER TABLE coding_sequences ADD UNIQUE (genbank_cds_id);
ALTER TABLE coding_sequences ADD UNIQUE (cds_code);
"""
)
conn.commit()

conn.close()

