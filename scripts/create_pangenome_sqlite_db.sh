#!/bin/bash
cd ${database}

## install relevant packages
sudo apt-get install sqlite

dbname=$1

### if no set db name env variable, INTERACTIVE prompt asks for database name and password
while [ -z $dbname ] ; do
  read -p 'Set PostgreSQL database name: ' dbname ; echo ''
done

## create database and load database schema
sqlite3 $dbname < $dbscripts/panterodb_initiate.sql

## populate database
# use tail to truncate the header
tail -n +2 ${allcomplete}_metadata_curated.tab > genome_assemblies.tab
cut -f1,3,4,5,6 ${allcomplete}_allreplicons_info.tab | tail -n +2 > genome_replicons.tab
cut -f1,2,3,4,5,6,8 ${allcomplete}_allproteins_info.tab | tail -n +2 > genome_coding_sequences.tab
cut -f2,3 $protali/full_families_info-noORFans.tab > genome_gene_families.tab
cut -f1,7 ${allcomplete}_allproteins_info.tab | tail -n +2 | grep -vP "^\t" | sort -u > genome_protein_products.tab 

sqlite3 $dbname << EOF
-- populate assemblies table
COPY assemblies FROM '$PWD/genome_assemblies.tab' NULL '';

-- populate replicons table
COPY replicons (assembly_id, genomic_accession, replicon_name, replicon_type, replicon_size) FROM '$PWD/genome_replicons.tab' NULL '';

-- populate protein table
begin;
create temp table protein_products (genbank_nr_protein_id VARCHAR(20), product TEXT) on commit drop;
create temp table protein_fams (genbank_nr_protein_id VARCHAR(20), protein_family_id CHAR(13)) on commit drop;
copy protein_products (genbank_nr_protein_id, product) from '$PWD/genome_protein_products.tab' NULL '';
copy protein_fams (protein_family_id, genbank_nr_protein_id) from '$protfamseqs.tab' NULL '';
insert into proteins (genbank_nr_protein_id, product, protein_family_id) (select genbank_nr_protein_id, min(product), protein_family_id
 from protein_products
 inner join protein_fams using (genbank_nr_protein_id)
 group by genbank_nr_protein_id,protein_family_id
);
commit;

-- populate coding_sequences table
begin;
create temp table codingsequences (
  nr_protein_id CHAR(15),
  genomic_accession CHAR(14) NOT NULL,
  locus_tag VARCHAR(200),
  cds_begin INTEGER NOT NULL,
  cds_end INTEGER NOT NULL,
  strand CHAR(1) NOT NULL,
  genbank_cds_id VARCHAR(50) NOT NULL
) on commit drop;
copy codingsequences from '$PWD/genome_coding_sequences.tab' null '';
create temp table cdsfam (
  genbank_cds_id VARCHAR(50) NOT NULL,
  gene_family_id CHAR(13)
) on commit drop;
copy cdsfam from '$PWD/genome_gene_families.tab';
create index gbcdsid_1 on codingsequences (genbank_cds_id);
create index gbcdsid_2 on cdsfam (genbank_cds_id);
insert into coding_sequences (genbank_cds_id, genomic_accession, locus_tag, cds_begin, cds_end, strand, genbank_nr_protein_id, gene_family_id) (
 select genbank_cds_id, genomic_accession, locus_tag, cds_begin, cds_end, strand, nr_protein_id, gene_family_id
  from codingsequences
  left join cdsfam using (genbank_cds_id)
);
commit;

-- populate the protein families
insert into nr_protein_families (select distinct protein_family_id, false from proteins);

-- add the bit mark for the singleton nr protein family
update nr_protein_families set is_singleton=true where protein_family_id='$protorfanclust';

-- allocated the '$cdsorfanclust' value to gene_family_id field for those CDSs with a parent protein but no gene family affiliation
update coding_sequences set gene_family_id='$cdsorfanclust' where gene_family_id is null and genbank_nr_protein_id is not null;

-- populate the gene families deriving from the protein families
insert into gene_families (select replace(protein_family_id, 'P', 'C') as gene_family_id, false, protein_family_id from nr_protein_families);

-- add the gene families which have no natural parent nr protein family, i.e. those derived from singleton nr proteins, but with multiple CDS members
insert into gene_families (select distinct gene_family_id, false, '$protorfanclust' from coding_sequences left join gene_families using (gene_family_id) where protein_family_id is null and gene_family_id is not null);

-- add the bit mark for the ORFan gene family
update gene_families set is_orfan=true where gene_family_id='$cdsorfanclust';

alter table gene_families add primary key (gene_family_id);

EOF

# load UniProt taxon codes for CDS name shortening
wget http://www.uniprot.org/docs/speclist
python << EOF
import re
import sqlite3, sys
conn = sqlite3.connect(database='$dbname')
cur = conn.cursor()
fout = open("uniprotcode_taxid.tab", 'w')
CODepat = re.compile("([A-Z0-9]{3,5}) +[ABEVO] +([0-9]{1,7}): ")
fspeclist = open('speclist', 'r')
for line in fspeclist:
  if line[0]==' ': continue
  CODEmatch = CODepat.match(line)
  if CODEmatch:
    code, taxid = CODEmatch.groups()
    fout.write("%s\t%s\n"%(code, taxid))

fout.close()
fspeclist.close()
cur.execute("COPY uniptrotcode2taxid FROM '$PWD/uniprotcode_taxid.tab' NULL '';")

# generate unique code for each genome assembly from the Uniprot code + enough digits
cur.execute("select assembly_id, uniptrotcode2taxid.code, species from assemblies left join uniptrotcode2taxid using (taxid);")
lasscode = cur.fetchall()
try:
  cur.execute("alter table assemblies add column code varchar(10);")
except sqlite3.ProgrammingError, e:
  conn.rollback()
  sys.stderr.write( 'Warning:' + str(e) + "reset assemblies.code values to NULL before update\n" )
  cur.execute("update assemblies set code=NULL;")

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
conn.close()
EOF

# add the cds_code column to coding_sequences table, generate unique cds codes based on the species code and genbank CDS id and create indexes
sqlite3 $dbname << EOF
create table coding_sequences2 (like coding_sequences);
alter table coding_sequences2 add column cds_code varchar(20);
insert into coding_sequences2 (
  select coding_sequences.*, code || '_' || regexp_replace(genbank_cds_id, '.+_([0-9]+)$', '\1') as cds_code
    from coding_sequences
    inner join replicons using (genomic_accession)
    inner join assemblies using (assembly_id)
  )
;
drop table coding_sequences;
alter table coding_sequences2 rename to coding_sequences;
-- recreate indexes and views
create index on coding_sequences (genbank_cds_id);
create index on coding_sequences (cds_code);
create index on coding_sequences (gene_family_id);
create index on coding_sequences (genomic_accession);
create index on coding_sequences (genomic_accession, cds_begin);
create index on coding_sequences (genbank_nr_protein_id);
create view gene_family_sizes as select gene_family_id, count(*) as size from coding_sequences group by gene_family_id;
alter table coding_sequences add UNIQUE (genbank_cds_id);
alter table coding_sequences add UNIQUE (cds_code);
EOF
