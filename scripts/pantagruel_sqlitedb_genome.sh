#!/bin/bash
thisscript=${0}
database=${1}
dbname=${2}
metadata=${3}
assemblyinfo=${4}
protali=${5}
protfamseqtab=${6}
protorfanclust=${7}
cdsorfanclust=${8}
usergenomeinfo=${9}
usergenomefinalassdir=${10}
gp2ass=${11}

echo "currently set variables:"
echo "database=${1} dbname=${2} metadata=${3} assemblyinfo=${4} protali=${5} protfamseqtab=${6} protorfanclust=${7} cdsorfanclust=${8} usergenomeinfo=$9 usergenomefinalassdir=${10} gp2ass=${11}"

if [ -z $cdsorfanclust ] ; then
 echo "Error: incomplete argument list. Usage:"
 echo "${thisscript} database dbname metadata assemblyinfo protali protfamseqtab protorfanclust cdsorfanclust [usergenomeinfo] [usergenomefinalassdir] [gp2ass]" exit 1
fi


initiatescript=${thisscript%.*}_initiate.sql
populatescript=${thisscript%.*}_populate.py

cd ${database}

# fetch UniProt taxon codes
if [ ! -e speclist ] ; then
 wget --progress=dot:mega http://www.uniprot.org/docs/speclist
fi


## create database and load database schema
sqlite3 ${dbname} < $initiatescript

## generate input data tables
step1="generating input files for database tables"
echo ${step1}
# create function for searching field positions
getnumfield (){
 local table=$1
 shift
 for f in ${@} ; do
  k=$(head -n 1 $table |  tr '\t' '\n' | grep -n $f)
  echo $k
 done | cut -d':' -f 1 | xargs | tr ' ' ','
}
# verify occurence of 'nr_protein_id' or 'protein_id' field in 'allproteins_info.tab'
protidfield=$(head -n 1 ${assemblyinfo}/allproteins_info.tab | grep -o 'nr_protein_id')
if [ -z "${protidfield}" ] ; then 
 protidfield=$(head -n 1 ${assemblyinfo}/allproteins_info.tab | grep -o 'protein_id')
 if [ -z "${protidfield}" ] ; then 
  echo "Error, file '${assemblyinfo}/allproteins_info.tab' if missing a field headed 'nr_protein_id' or 'protein_id'"
  exit 1
 fi
fi
protidsedfilter='s/[^\t]*protein_id[^\t]*/nr_protein_id/'
# assembly table
cat ${metadata}/metadata_curated.tab > genome_assemblies.tab
# replicon table
replifields="assembly_id genomic_accession replicon_name replicon_type replicon_size"
replifieldnums=$(getnumfield ${assemblyinfo}/allreplicons_info.tab $replifields)
cut -f $replifieldnums ${assemblyinfo}/allreplicons_info.tab > genome_replicons.tab
# gene_family table
cut -f 2,3 ${protali}/full_families_info-noORFans.tab > genome_gene_families.tab
# protein_family table
awk -v"OFS=\t" '{print $2 ,$1}' ${protfamseqs}.tab ${protfamseqtab} > genome_protein_families.tab
# coding_sequences table
cdsfields="$protidfield genomic_accession locus_tag cds_begin cds_end cds_strand genbank_cds_id"
cdsfieldnums=$(getnumfield ${assemblyinfo}/allproteins_info.tab $cdsfields)
head -n 1  ${assemblyinfo}/allproteins_info.tab | cut -f ${cdsfieldnums} | sed -e ${protidsedfilter} > genome_coding_sequences.tab && \
tail -n +2 ${assemblyinfo}/allproteins_info.tab | cut -f ${cdsfieldnums} >> genome_coding_sequences.tab
# proteins table
protfields="$protidfield product"
protfieldnums=$(getnumfield ${assemblyinfo}/allproteins_info.tab $protfields)
head -n 1  ${assemblyinfo}/allproteins_info.tab | cut -f ${protfieldnums} | sed -e ${protidsedfilter} > genome_protein_products.tab && \
tail -n +2 ${assemblyinfo}/allproteins_info.tab | cut -f ${protfieldnums} | grep -vP "^\t" | sort -u | sed -e "s/%2C/,/g" | sed -e "s/%3B/;/g" >> genome_protein_products.tab

ls -l genome_replicons.tab genome_gene_families.tab  genome_protein_families.tab genome_coding_sequences.tab genome_protein_products.tab
checkexec "failed ${step1}" "succesfully ${step1/ing/ed}"

## populate database
step2="populating database with information on genome assemblies, CDS/protein annotation and gene families"
echo ${step2}
if [ -z  "${customassemb}" ] ; then
  python2.7 "${populatescript}" -d "${dbname}" -o "${protorfanclust}" -O "${cdsorfanclust}" -s "${database}/speclist" -i "${gp2ass}"
else
  python2.7 "${populatescript}" -d "${dbname}" -o "${protorfanclust}" -O "${cdsorfanclust}" -s "${database}/speclist" -i "${gp2ass}" --user.genome.info "${usergenomeinfo}" --user.genome.ass.dir "${usergenomefinalassdir}"
fi
checkexec "failed ${step2}" "succesfully ${step2/ing/ed}"
cd -
