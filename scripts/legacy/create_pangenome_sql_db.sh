#!/bin/bash
cd ${database}

## install relevant packages
sudo apt-get install postgresql postgresql-contrib python-psycopg2
# install graphic interfaces for local and distant (web-based) admin
sudo apt-get install pgadmin phppgadmin

dbname=$1

### INTERACTIVE prompts ask for database name and password
while [ -z $dbname ] ; do
  read -p 'Set PostgreSQL database name: ' dbname ; echo ''
done

## check existence of a superuser with same name as current linux user
if [ -z "$(sudo -u postgres psql -c 'SELECT usename FROM pg_user;' | grep $USER)" ] ; then
# if not create it
while [ -z $passvar ] ; do
  read -sp "Set PostgreSQL password for user '$USER': " passvar ; echo ''
  read -sp "Confirm password: " passvar2 ; echo ''
  if [[ "$passvar" != "$passvar2" ]] ; then
    echo "Entered passwords are empty or do not match; please enter password again"
    unset passvar passvar2
  fi
done
unset passvar2
sudo -u postgres createuser --superuser --echo -P << EOF
$passvar
$passvar
EOF
# create pgpass file
echo "localhost:5432:$dbname:$USER:$passvar" >> ~/.pgpass
chmod 0600 ~/.pgpass
unset passvar 
fi
## create database
sudo -u postgres createdb $dbname
# add extension for use of pgadmin and change database ownership
sudo -u postgres psql << EOF
CREATE EXTENSION adminpack;
ALTER DATABASE "$dbname" OWNER TO $USER;
\q
EOF

## load database schema
psql -h localhost -d "$dbname" < $dbscripts/panterodb_initiate.sql

## populate database
# use tail to truncate the header
#~ # insert an emptyfield at the end of 'genome_assemblies.tab' to account for 'code' column, to be filled later.
#~ tail -n +2 ${allcomplete}_metadata_curated.tab | sed -e 's/$/\t/g' > genome_assemblies.tab
tail -n +2 ${allcomplete}_metadata_curated.tab > genome_assemblies.tab
cut -f1,3,4,5,6 ${allcomplete}_allreplicons_info.tab | tail -n +2 > genome_replicons.tab
cut -f1,2,3,4,5,6,8 ${allcomplete}_allproteins_info.tab | tail -n +2 > genome_coding_sequences.tab
cut -f2,3 $protali/full_families_info-noORFans.tab > genome_gene_families.tab
cut -f1,7 ${allcomplete}_allproteins_info.tab | tail -n +2 | grep -vP "^\t" | sort -u > genome_protein_products.tab 

# populate assemblies table
psql -h localhost -d "$dbname" -c "COPY genome.assemblies FROM '$PWD/genome_assemblies.tab' NULL '';"
#~ COPY 880
# populate replicons table
psql -h localhost -d "$dbname" -c "COPY genome.replicons (assembly_id, genomic_accession, replicon_name, replicon_type, replicon_size) FROM '$PWD/genome_replicons.tab' NULL '';"
#~ COPY 2361
# populate protein table
psql -h localhost -d "$dbname" << EOF
begin;
create temp table protein_products (genbank_nr_protein_id VARCHAR(20), product TEXT) on commit drop;
create temp table protein_fams (genbank_nr_protein_id VARCHAR(20), protein_family_id CHAR(13)) on commit drop;
copy protein_products (genbank_nr_protein_id, product) from '$PWD/genome_protein_products.tab' NULL '';
copy protein_fams (protein_family_id, genbank_nr_protein_id) from '$protfamseqs.tab' NULL '';
insert into genome.proteins (genbank_nr_protein_id, product, protein_family_id) (select genbank_nr_protein_id, min(product), protein_family_id
 from protein_products
 inner join protein_fams using (genbank_nr_protein_id)
 group by genbank_nr_protein_id,protein_family_id
);
commit;
EOF
#~ BEGIN
#~ CREATE TABLE
#~ CREATE TABLE
#~ COPY 849727
#~ COPY 834474
#~ INSERT 0 834474
#~ COMMIT

# populate coding_sequences table
psql -h localhost -d "$dbname" << EOF
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
insert into genome.coding_sequences (genbank_cds_id, genomic_accession, locus_tag, cds_begin, cds_end, strand, genbank_nr_protein_id, gene_family_id) (
 select genbank_cds_id, genomic_accession, locus_tag, cds_begin, cds_end, strand, nr_protein_id, gene_family_id
  from codingsequences
  left join cdsfam using (genbank_cds_id)
);
commit;
EOF
#~ BEGIN
#~ CREATE TABLE
#~ COPY 4505870
#~ CREATE TABLE
#~ COPY 4301813
#~ CREATE INDEX
#~ CREATE INDEX
#~ INSERT 0 4505870
#~ COMMIT

# populate the protein families
psql -h localhost -d "$dbname" -c "insert into genome.nr_protein_families (select distinct protein_family_id, false from genome.proteins);"
#~ INSERT 0 22743
# add the bit mark for the singleton nr protein family
psql -h localhost -d "$dbname" -c "update genome.nr_protein_families set is_singleton=true where protein_family_id='$protorfanclust';"
#~ UPDATE 1
# allocated the '$cdsorfanclust' value to gene_family_id field for those CDSs with a parent protein but no gene family affiliation
psql -h localhost -d "$dbname" -c "update genome.coding_sequences set gene_family_id='$cdsorfanclust' where gene_family_id is null and genbank_nr_protein_id is not null;"
#~ UPDATE 14771
# populate the gene families deriving from the protein families
psql -h localhost -d "$dbname" -c "insert into genome.gene_families (select replace(protein_family_id, 'P', 'C') as gene_family_id, false, protein_family_id from genome.nr_protein_families);"
#~ INSERT 0 22743
# add the gene families which have no natural parent nr protein family, i.e. those derived from singleton nr proteins, but with multiple CDS members
psql -h localhost -d "$dbname" -c "insert into genome.gene_families (select distinct gene_family_id, false, '$protorfanclust' from genome.coding_sequences left join genome.gene_families using (gene_family_id) where protein_family_id is null and gene_family_id is not null);"
#~ INSERT 0 4424
# add the bit mark for the ORFan gene family
psql -h localhost -d "$dbname" -c "update genome.gene_families set is_orfan=true where gene_family_id='$cdsorfanclust';"
#~ UPDATE 1
psql -h localhost -d "$dbname" -c "alter table genome.gene_families add primary key (gene_family_id);"
#~ ALTER TABLE

# load UniProt taxon codes for CDS name shortening
wget http://www.uniprot.org/docs/speclist
python << EOF
import re
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
EOF
psql -h localhost -d "$dbname" -c "COPY taxonomy.uniptrotcode2taxid FROM '$PWD/uniprotcode_taxid.tab' NULL '';"

# generate unique code for each genome assembly from the Uniprot code + enough digits
python << EOF
import psycopg2, sys
conn = psycopg2.connect(dbname='$dbname', user='ubuntu', host='localhost')
cur = conn.cursor()
cur.execute("select assembly_id, uniptrotcode2taxid.code, species from genome.assemblies left join taxonomy.uniptrotcode2taxid using (taxid);")
lasscode = cur.fetchall()
try:
  cur.execute("alter table genome.assemblies add column code varchar(10);")
except psycopg2.ProgrammingError, e:
  conn.rollback()
  sys.stderr.write( 'Warning:' + str(e) + "reset assemblies.code values to NULL bfore update\n" )
  cur.execute("update genome.assemblies set code=NULL;")

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
  cur.execute("update genome.assemblies set code=%s where assembly_id=%s;", (c+str(dcodesn[c]), ass.strip(' ')))
  print c+str(dcodesn[c]), ass
conn.commit()
EOF

# add the cds_code column to genome.coding_sequences table, generate unique cds codes based on the species code and genbank CDS id and create indexes
psql -h localhost -d "$dbname" << EOF
set search_path = genome, pg_catalog;
create table coding_sequences2 (like genome.coding_sequences);
alter table genome.coding_sequences2 add column cds_code varchar(20);
insert into genome.coding_sequences2 (
  select coding_sequences.*, code || '_' || regexp_replace(genbank_cds_id, '.+_([0-9]+)$', '\1') as cds_code
    from genome.coding_sequences
    inner join genome.replicons using (genomic_accession)
    inner join genome.assemblies using (assembly_id)
  )
;
drop table genome.coding_sequences;
alter table genome.coding_sequences2 rename to coding_sequences;
-- recreate indexes and views
create index on genome.coding_sequences (genbank_cds_id);
create index on genome.coding_sequences (cds_code);
create index on genome.coding_sequences (gene_family_id);
create index on genome.coding_sequences (genomic_accession);
create index on genome.coding_sequences (genomic_accession, cds_begin);
create index on genome.coding_sequences (genbank_nr_protein_id);
create view gene_family_sizes as select gene_family_id, count(*) as size from genome.coding_sequences group by gene_family_id;
alter table genome.coding_sequences add UNIQUE (genbank_cds_id);
alter table genome.coding_sequences add UNIQUE (cds_code);
EOF
