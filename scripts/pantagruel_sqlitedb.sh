#!/bin/bash
thisscript=$0
database=$1
dbname=$2
metadata=$3
assemblyinfo=$4
protali=$5
protfamseqtab=$6
protorfanclust=$7
cdsorfanclust=$8

initiatescript=${thisscript%.*}_initiate.sql
populatescript=${thisscript%.*}_populate.py

cd ${database}

# fetch UniProt taxon codes
if [ ! -e speclist ] ; then
 wget http://www.uniprot.org/docs/speclist
fi
### if no set db name env variable, INTERACTIVE prompt asks for database name and password
while [ -z $dbname ] ; do
  read -p 'Set PostgreSQL database name: ' dbname ; echo ''
done

## create database and load database schema
sqlite3 $dbname < $initiatescript

## populate database
# test presence of 'nr_protein_id' field in allproteins_info.tab (should be 9th field)
protidfield=$(head -n 1 ${assemblyinfo}/allproteins_info.tab |  tr '\t' '\n' | grep -n 'nr_protein_id' | cut -d':' -f1)
if [ -z $protidfield ] ; then 
 # otherwise select 'protein_id' field (should be 1st field)
 protidfield=$(head -n 1 ${assemblyinfo}/allproteins_info.tab |  tr '\t' '\n' | grep -n 'protein_id' | cut -d':' -f1)
fi
# use tail to truncate the header
cat ${metadata}/metadata_curated.tab > genome_assemblies.tab
cut -f1,3,4,5,6 ${assemblyinfo}/allreplicons_info.tab > genome_replicons.tab
cut -f${protidfield},2,3,4,5,6,8 ${assemblyinfo}/allproteins_info.tab > genome_coding_sequences.tab
cut -f2,3 ${protali}/full_families_info-noORFans.tab > genome_gene_families.tab
cut -f${protidfield},7 ${assemblyinfo}/allproteins_info.tab | grep -vP "^\t" | sort -u > genome_protein_products.tab 

# populate database
python $populatescript $dbname ${protfamseqtab} ${protorfanclust} ${cdsorfanclust}

cd -
