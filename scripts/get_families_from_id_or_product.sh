#!/bin/bash
genefamlist=${1}
allfamlist=${2}

for fam in $(cut -f1 ${genefamlist} | sed -e 's/ /=^+^=/g') ; do
  if [ ! -z "$(echo ${fam} | grep ${famprefix})" ] ; then
    # a gene family id
    grep ${fam} ${allfamlist}
  else
    # anything else, assume its a gene product description
	prod="$(echo ${fam} | sed -e 's/=^+^=/ /g' | sed -e "s/'/\'/g")"
	sqlite3 ${sqldb} "select distinct gene_family_id from coding_sequences inner join proteins using (nr_protein_id) where product like '%${fam}%';"
  fi
done