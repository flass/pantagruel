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
sqlite3 $dbname < $dbscripts/pantagruel_sqlitedb_initiate.sql

## populate database
# use tail to truncate the header
tail -n +2 ${allcomplete}_metadata_curated.tab > genome_assemblies.tab
cut -f1,3,4,5,6 ${allcomplete}_allreplicons_info.tab | tail -n +2 > genome_replicons.tab
cut -f1,2,3,4,5,6,8 ${allcomplete}_allproteins_info.tab | tail -n +2 > genome_coding_sequences.tab
cut -f2,3 $protali/full_families_info-noORFans.tab > genome_gene_families.tab
cut -f1,7 ${allcomplete}_allproteins_info.tab | tail -n +2 | grep -vP "^\t" | sort -u > genome_protein_products.tab 

python $dbscripts/pantagruel_sqlitedb_populate.py $dbname

# load UniProt taxon codes for CDS name shortening
wget http://www.uniprot.org/docs/speclist
python << EOF
import re
import sqlite3, sys
conn = sqlite3.connect(database='$dbname')
cur = conn.cursor()
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
