#!/bin/bash


inrefass=${1}
refgenus=${2}
if [ -z ${3} ] ; then
  tmpdir=${3}
else
  tmpdir=$PWD
fi
if [ -z ${4} ] ; then
  logdir=${4}
else
  logdir=$PWD
fi

prokkabin=$(which prokka)
prokkadir=$(dirname $(dirname $(readlink -f $prokkabin)))
prokkablastdb=${prokkadir}/db/genus

testpgb2fdb () {
  gbff=${1}
  testgbff=${gbff%*.gz}
  if [ "${testgbff}" != "${gbff}" ] ; then gunzip -k ${gbff} ; fi
  # test the functioning of prokka-genbank_to_fasta_db script
  pgb2fdb="prokka-genbank_to_fasta_db"
  ${pgb2fdb} ${testgbff} &> /dev/null
  if [ ${?} -gt 0 ] ; then
    # try using the native perl
    pgb2fdb="/usr/bin/perl $prokkadir/bin/prokka-genbank_to_fasta_db"
    ${pgb2fdb} ${testgbff} > /dev/null
    if [ ${?} -gt 0 ] ; then
    echo "Error: cannot make prokka-genbank_to_fasta_db to work (likely a problem of Perl lib dependencies) ; exit now" >&2
    exit 1
    fi
  fi
  if [ "${testgbff}" != "${tgbff}" ] ; then rm -f ${testgbff} ; fi
  echo ${pgb2fdb}
}
export -f testpgb2fdb

# function for extraction
gzpgb2fdb (){
  pgb2fdb=${1}
  nfgffgz=${2}
  nfout=${3}
  nflog=${4}
  gunzip -k ${nfgffgz}
  nfgff=${nfgffgz%*.gz}
  echo "${pgb2fdb} ${nfgff}  ..."
  eval "${pgb2fdb} ${nfgff} 1> ${nfout}.tmp 2>> ${nflog}"
  rm -f ${nfgff}
  if [ -s ${nfout}.tmp ] ; then
    echo "${pgb2fdb} ${nfgff}  ... done"
    echo "output:"
    ls -l ${nfout}.tmp
    cat ${nfout}.tmp >> ${nfout} && rm -f ${nfout}.tmp
	return 0
  else
#    return 1
    # for the moment produces a warning only as a bug in prokka
	# results in recurent failure to produce CDSs from GenBank flat files
	# cf. https://github.com/tseemann/prokka/issues/ #400 and #393
    echo "WARNING: produced no CDS for file $nfgff" >&2
    echo "${pgb2fdb} ${nfgff}  ... failed"
  fi
}
export -f gzpgb2fdb


if [ -e ${prokkablastdb}/${refgenus} ] ; then
  datetime=$(date +%Y-%m-%d_%H-%M-%S)
  echo "Warning: found a previous Prokka reference database for the genus '${refgenus}; save it to '${prokkablastdb}/${refgenus}.backup_${datetime}'"
  mv -f ${prokkablastdb}/${refgenus} ${prokkablastdb}/${refgenus}.backup_${datetime}
fi
# add all the (representative) proteins in the dataset to the custom reference prot database for Prokka to search for similarities
ls ${inrefass}/*/*_genomic.gbff.gz > ${inrefass}_genomic_gbffgz_list
pgb2fdb=$(testpgb2fdb `head -n1 ${inrefass}_genomic_gbffgz_list`)
rm -f ${tmpdir}/${refgenus}.faa.*  ${logdir}/prokka-genbank_to_fasta_db.log.*
if [ -s ${inrefass}_genomic_gbffgz_list ] ; then
  # extract protein sequences
  parallel --halt now,fail=1 -a ${inrefass}_genomic_gbffgz_list gzpgb2fdb "${pgb2fdb} {} ${tmpdir}/${refgenus}.faa.{%} ${logdir}/prokka-genbank_to_fasta_db.log.{%}"
  if [ ${?} -gt 0 ] ; then
    >&2 echo "WARNING: could not build custom genus-specific BLAST db from files ${inrefass}/*/*_genomic.gbff ; Prokka will have to use default databases."
  else
    cat ${tmpdir}/${refgenus}.faa.* > ${tmpdir}/${refgenus}.faa && rm -f ${tmpdir}/${refgenus}.faa.*
    # cluster similar sequences
    cdhit -i ${tmpdir}/${refgenus}.faa -o ${tmpdir}/${refgenus}_representative.faa -T 0 -M 0 -G 1 -s 0.8 -c 0.9 &> ${logdir}/cdhit.log
    rm -fv ${tmpdir}/${refgenus}.faa ${tmpdir}/${refgenus}_representative.faa.clstr
    # replace database name for genus detection by Prokka 
    cp -p ${tmpdir}/${refgenus}_representative.faa ${prokkablastdb}/${refgenus}
    cd ${prokkablastdb}/
    makeblastdb -dbtype prot -in ${refgenus}
    cd -
  fi          
fi
