#!/usr/bin/python
import re
import os, sys
import sqlite3

def getTableInfos(table, cursor, ommitserial=False):
	cursor.execute("PRAGMA table_info(%s);"%table)
	colinfos = cursor.fetchall()
	if ommitserial:
		 colinfos = [colinfo for colinfo in colinfos if colinfo['type'].upper()!='SERIAL']
	return colinfos 

def replaceValuesAsNull(table, cursor, nullval='', ommitcols=[], tableinfo=None):
	if not tableinfo: colinfos = getTableInfos(table, cursor)
	else: colinfos = tableinfo
	for colinfo in colinfos:
		colname = colinfo['name']
		if (colname not in ommitcols):
			cursor.execute("UPDATE %s SET %s=NULL WHERE %s='%s';"%(table, colname, colname, nullval))

def loadAndCurateTable(table, nfin, cursor, header=True, insertcolumns=(), sep='\t', doNotReplaceWithNull=[], **kw):
	ftabin = open(nfin, 'r')
	colinfos = getTableInfos(table, cursor, ommitserial=True)
	insertcols=[]
	if header:
		insertcols = ftabin.readline().rstrip('\n').split(sep)
	elif insertcolumns:
		if type(insertcolumns)==str: insertcols = [incol.strip(' ') for incol in insertcolumns.split(',')]
		else: insertcols = insertcolumns
	else:
		insertcols = [colinfo['name'] for colinfo in colinfos]
	coldef = '('+', '.join(insertcols)+')'
	print table, coldef
	cursor.executemany("INSERT INTO %s %s VALUES (%s);"%(table, coldef, ','.join(['?']*len(insertcols))), (tuple(line.rstrip('\n').split(sep)) for line in ftabin))
	replaceValuesAsNull(table, cursor, tableinfo=colinfos, ommitcols=doNotReplaceWithNull, **kw)
	ftabin.close()
	
def createAndLoadTable(table, tabledef, nfin, cursor, temp=False, enddrop=False, **kw):
	tmp = 'TEMP' if temp else ''
	if tabledef.upper().startswith('LIKE ') and 'header' in kw:
		with open(nfin, 'r') as ftabin:
			insertcols = ftabin.readline().rstrip('\n').split(sep)
		reftable = tabledef.split('LIKE ', 1)[1]
		refcolinfos = getTableInfos(reftable, cursor, ommitserial=True)
		tdeflines = []
		for incol in insertcols:
			for refcolinfo in refcolinfos:
				if refcolinfo['name']==incol:
					tdefline = [refcolinfo['name'], 
					            refcolinfo['type'], 
					            'NOT NULL' if bool(int(refcolinfo['nn'])) else 'NULL', 
					            'DEFAULT'+refcolinfo['dflt_value'] if refcolinfo['dflt_value'] else '', 
					            'PRIMARY KEY'  if bool(int(refcolinfo['pk'])) else '']
					tdeflines.append(' '.join(tdefline))
					break
			else:
				raise ValueError, "no collumn named '%s' in reference table '%s'"%(incol, reftable)
		tdef = '('+', '.join(tdeflines)+')'
	else:
		tdef = tabledef
	cursor.execute("CREATE %s TABLE %s %s;"%(tmp, table, tdef))
	loadAndCurateTable(table, nfin, cursor, **kw)
	if enddrop: cursor.execute("DROP TABLE %s;"%table)


cdsnumpat = re.compile('.+_([0-9]+)$')
def make_cds_code(code, genbank_cds_id):
	return "%s_%s"%(code, cdsnumpat.search(genbank_cds_id).group(1))


dbname = sys.argv[1] # os.environ['sqldbname']
protfamseqtab = sys.argv[2] # os.environ['protfamseqs']+'.tab'
protorfanclust = sys.argv[3] # os.environ['protorfanclust']
cdsorfanclust = sys.argv[4] # os.environ['cdsorfanclust']
nfspeclist = sys.argv[5] # 'speclist'
if len(sys.argv) > 6:
	nfusergenomeinfo = sys.argv[6]
else:
	nfusergenomeinfo = None

conn = sqlite3.connect(database=dbname)
conn.create_function("make_cds_code", 2, make_cds_code)
conn.row_factory = sqlite3.Row
cur = conn.cursor()

# populate assemblies table
loadAndCurateTable('assemblies', 'genome_assemblies.tab', cur, header=True, doNotReplaceWithNull=['assembly_id', 'assembly_name', 'taxid'])

# create indexes on assemblies table
cur.executescript("""
CREATE INDEX assemblies_assembly_id_key ON assemblies (assembly_id);
CREATE INDEX assemblies_species_key ON assemblies (species);
CREATE INDEX assemblies_taxid_key ON assemblies (taxid);
""")
conn.commit()

# populate replicons table
#~ loadAndCurateTable('replicons', 'genome_replicons.tab', cur, insertcolumns="(assembly_id, genomic_accession, replicon_name, replicon_type, replicon_size)")
loadAndCurateTable('replicons', 'genome_replicons.tab', cur, header=True)
conn.commit()

# populate protein table
protprodtabledef = """(nr_protein_id VARCHAR(20), product TEXT)"""
protfamtabledef = """(nr_protein_id VARCHAR(20), protein_family_id CHAR(13))"""
createAndLoadTable('protein_products', protprodtabledef, 'genome_protein_products.tab', cur)
createAndLoadTable('protein_fams', protfamtabledef, protfamseqtab, cur, header=False)
cur.executescript("""
INSERT INTO proteins (nr_protein_id, product, protein_family_id)
 SELECT nr_protein_id, min(product), protein_family_id
  FROM protein_products
  INNER JOIN protein_fams USING (nr_protein_id)
  GROUP BY nr_protein_id,protein_family_id
;
DROP TABLE protein_products;
DROP TABLE protein_fams;
""")
conn.commit()

# create indexes on proteins table
cur.executescript("""
CREATE INDEX proteins_nr_protein_id_key ON proteins (nr_protein_id);
CREATE INDEX proteins_protein_family_id_key ON proteins (protein_family_id);
CREATE INDEX proteins_product_key ON proteins (product);
""")
conn.commit()

# populate coding_sequences table
cdstabledef = """(
  nr_protein_id CHAR(15),
  genomic_accession CHAR(14) NOT NULL,
  locus_tag VARCHAR(200),
  cds_begin INTEGER NOT NULL,
  cds_end INTEGER NOT NULL,
  cds_strand CHAR(1) NOT NULL,
  genbank_cds_id VARCHAR(50) NOT NULL
)"""
cdsfamtabledef = """(
  genbank_cds_id VARCHAR(50) NOT NULL,
  gene_family_id CHAR(13)
)"""
#~ createAndLoadTable('codingsequences', 'LIKE coding_sequences', 'genome_coding_sequences.tab', cur, header=True)
#~ createAndLoadTable('cdsfam', 'LIKE coding_sequences', 'genome_gene_families.tab', cur, header=True)
createAndLoadTable('codingsequences', cdstabledef, 'genome_coding_sequences.tab', cur, header=True)
createAndLoadTable('cdsfam', cdsfamtabledef, 'genome_gene_families.tab', cur, header=False)
cur.executescript("""
CREATE INDEX gbcdsid_1 ON codingsequences (genbank_cds_id);
CREATE INDEX gbcdsid_2 ON cdsfam (genbank_cds_id);
INSERT INTO coding_sequences (genbank_cds_id, genomic_accession, locus_tag, cds_begin, cds_end, cds_strand, nr_protein_id, gene_family_id) 
 SELECT genbank_cds_id, genomic_accession, locus_tag, cds_begin, cds_end, cds_strand, nr_protein_id, gene_family_id
  FROM codingsequences
  LEFT JOIN cdsfam USING (genbank_cds_id);
DROP TABLE codingsequences;
DROP TABLE cdsfam;
""")
conn.commit()

# populate the protein families
cur.execute("INSERT INTO nr_protein_families (protein_family_id) SELECT DISTINCT protein_family_id FROM proteins;")
# add the bit mark for the singleton nr protein family
cur.execute("UPDATE nr_protein_families SET is_singleton=1 WHERE protein_family_id=?;", (protorfanclust,))
# allocated the '$cdsorfanclust' value to gene_family_id field for those CDSs with a parent protein but no gene family affiliation
cur.execute("UPDATE coding_sequences SET gene_family_id=? " \
           +"WHERE gene_family_id IS NULL AND nr_protein_id IS NOT NULL;", (cdsorfanclust,))
# populate the gene families deriving from the protein families
cur.execute("INSERT INTO gene_families " \
           +"SELECT replace(protein_family_id, 'P', 'C') AS gene_family_id, 0, protein_family_id " \
           +"FROM nr_protein_families;")
# add the gene families which have no natural parent nr protein family, i.e. those derived from singleton nr proteins, but with multiple CDS members
cur.execute("INSERT INTO gene_families SELECT distinct gene_family_id, 0, ? " \
           +"FROM coding_sequences LEFT JOIN gene_families USING (gene_family_id) " \
           +"WHERE protein_family_id IS NULL AND gene_family_id IS NOT NULL;", (protorfanclust,))
# add the bit mark for the ORFan gene family
cur.execute("UPDATE gene_families SET is_orfan=1 WHERE gene_family_id=?;", (cdsorfanclust,))
conn.commit()

# create indexes on protein/gene family tables
cur.executescript("""
CREATE INDEX nrprotfams_protein_family_id_key ON nr_protein_families (protein_family_id);
CREATE INDEX genefams_gene_family_id_key ON gene_families (gene_family_id);
CREATE INDEX genefams_protein_family_id_key ON gene_families (protein_family_id);
""")
conn.commit()

# load UniProt taxon codes for CDS name shortening
fout = open("uniprotcode_taxid.tab", 'w')
codetaxids = []
CODepat = re.compile("([A-Z0-9]{3,5}) +[ABEVO] +([0-9]{1,7}): ")
with open(nfspeclist, 'r') as fspeclist:
	for line in fspeclist:
		if line[0] in [' ', '<']: continue # skip comments and HTML content
		CODEmatch = CODepat.match(line)
		if CODEmatch:
			code, taxid = CODEmatch.groups()
			fout.write("%s\t%s\n"%(code, taxid))
			codetaxids.append((code, taxid))

if nfusergenomeinfo:
	with open(nfusergenomeinfo, 'r') as fusergenomeinfo:
		uginfheader = fusergenomeinfo.readline().rstrip('\n').split('\t')
		for line in fusergenomeinfo:
			lsp = line.rstrip('\n').split('\t')
			duginfo = {field:lsp[f] for f, field in enumerate(uginfheader)}
			code, taxid = (duginfo['locus_tag_prefix'], duginfo['taxid'])
			fout.write("%s\t%s\n"%(code, taxid))
			codetaxids.append((code, taxid))

fout.close()
cur.executemany("INSERT INTO uniptrotcode2taxid VALUES (?, ?);", codetaxids)
conn.commit()

# generate unique code for each genome assembly from the Uniprot code + enough digits
cur.execute("SELECT assembly_id, uniptrotcode2taxid.code, species FROM assemblies LEFT JOIN uniptrotcode2taxid USING (taxid) ORDER BY uniptrotcode2taxid.code DESC;")
lasscode = cur.fetchall()
#~ try:
cur.execute("ALTER TABLE assemblies ADD COLUMN code VARCHAR(10) DEFAULT NULL;")
#~ except sqlite3.Error, e:
  #~ conn.rollback()
  #~ sys.stderr.write( 'Warning:' + str(e) + "reset assemblies.code values to NULL before update\n" )
  #~ cur.execute("UPDATE assemblies SET code=NULL;")

dcodesn = {}
dcodeass = {}
for ass, code, spe in lasscode:
	print ass, code, spe
	if code:
		c = str(code)
	else:
		spsp = str(spe).split()
		s = spsp[-2]
		p = spsp[-1]
		if '.' in p:
			# fairly common case of organism name being '[Candidatus] Genus sp.'
			c = s[:6].upper()
		else:
			# regular case  of organism name being '[Candidatus] Genus species'
			c = (s[:3]+p[:3]).upper()
	dcodesn[c] = dcodesn.setdefault(c, 0) + 1
	if dcodesn[c] == 1:
		dcodeass[c] = ass
		print c, ass
	else:
		dcodeass[c+'1'] = dcodeass[c]
		del dcodeass[c]
		dcodeass[c+str(dcodesn[c])] = ass
		print c+str(dcodesn[c]), ass

cur.executemany("UPDATE assemblies set code=? WHERE assembly_id=?;", dcodeass.iteritems())
conn.commit()

# add the cds_code column to coding_sequences table, generate unique cds codes based on the species code and genbank CDS id and create indexes
cur.executescript("""
CREATE TABLE coding_sequences2 AS SELECT * FROM coding_sequences WHERE 0;
ALTER TABLE coding_sequences2 ADD COLUMN cds_code varchar(20);
INSERT INTO coding_sequences2 
  SELECT coding_sequences.*, make_cds_code(code, genbank_cds_id) AS cds_code
    FROM coding_sequences
    INNER JOIN replicons USING (genomic_accession)
    INNER JOIN assemblies USING (assembly_id);
""")
conn.commit()

# drop original table and create indexes and views
cur.executescript("""
DROP TABLE coding_sequences;
ALTER TABLE coding_sequences2 RENAME TO coding_sequences;
CREATE INDEX cds_genbank_cds_id_key ON coding_sequences (genbank_cds_id);
CREATE INDEX cds_cds_code_key ON coding_sequences (cds_code);
CREATE INDEX cds_gene_family_id_key ON coding_sequences (gene_family_id);
CREATE INDEX cds_genomic_accession_key ON coding_sequences (genomic_accession);
CREATE INDEX cds_genomic_accession_cds_begin_key ON coding_sequences (genomic_accession, cds_begin);
CREATE INDEX cds_nr_protein_id_key ON coding_sequences (nr_protein_id);
CREATE VIEW gene_family_sizes AS SELECT gene_family_id, count(*) as size from coding_sequences group by gene_family_id;
""")
conn.commit()

conn.close()

