#!/bin/bash

database=$1
dbname=$2
assemblyinfo=$3
protali=$4
protfamseqtab=$5
protorfanclust=$6
cdsorfanclust=$7
initiatescript=$8
populatescript=$9

cd ${database}

# fetch UniProt taxon codes
wget http://www.uniprot.org/docs/speclist

### if no set db name env variable, INTERACTIVE prompt asks for database name and password
while [ -z $dbname ] ; do
  read -p 'Set PostgreSQL database name: ' dbname ; echo ''
done

## create database and load database schema
sqlite3 $dbname < $initiatescript

## populate database
# use tail to truncate the header
tail -n +2 ${assemblyinfo}/metadata_curated.tab > genome_assemblies.tab
cut -f1,3,4,5,6 ${assemblyinfo}/allreplicons_info.tab | tail -n +2 > genome_replicons.tab
cut -f1,2,3,4,5,6,8 ${assemblyinfo}/allproteins_info.tab | tail -n +2 > genome_coding_sequences.tab
cut -f2,3 ${protali}/full_families_info-noORFans.tab > genome_gene_families.tab
# test presence of 'nr_protein_id' field
protidfield=$(head -n 1 ${assemblyinfo}/allproteins_info.tab |  tr '\t' '\n' | | grep -n 'nr_protein_id' | cut -d':' -f1)
if [ -z $protidfield ] ; then 
 protidfield=$(head -n 1 ${assemblyinfo}/allproteins_info.tab |  tr '\t' '\n' | | grep -n 'protein_id' | cut -d':' -f1)
fi
cut -f${protidfield},7 ${assemblyinfo}/allproteins_info.tab | tail -n +2 | grep -vP "^\t" | sort -u > genome_protein_products.tab 

# populate database
python $populatescript $dbname ${protfamseqtab} ${protorfanclust} ${cdsorfanclust}
