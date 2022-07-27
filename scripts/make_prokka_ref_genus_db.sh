#!/bin/bash

echo "This script relies on GNU Parallel for multi-threading"
parallel --citation


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

## parse script arguments
inrefass=${1}
refgenus=${2}
if [ ! -z ${3} ] ; then
  tmpdir=${3}
else
  tmpdir=$PWD
fi
if [ ! -z ${4} ] ; then
  logdir=${4}
else
  logdir=$PWD
fi

# test existence of prokka command
prokkabin=$(which prokka)
if [ -z "${prokkabin}" ] ; then
  if [[ "${inrefass}" == 'check' ]] ; then
    echo "noprokka"
    exit 0
  else
    &>2 echo "Error: command 'prokka' was not found"
    exit 1
  fi
fi
prokkadir=$(dirname $(dirname $(readlink -f $prokkabin)))
if [ ! -z "$(env | grep ${PROKKA_DATA_DIR} | cut -d'=' -f2)" ] ; then
  prokkablastdb=${PROKKA_DATA_DIR}
else
  prokkablastdb=${prokkadir}/db/genus
fi

if [[ "${inrefass}" == 'check' ]] ; then
  # just check the pre-existence of the database and exit with code 0
  # print the path of the db file exists with non-null size ; nothing if the db file is absent or empty
  if [ -s "${prokkablastdb}/${refgenus}" ] ; then
    echo ${prokkablastdb}/${refgenus}
  elif [ ! -w "${prokkablastdb}/" ] ; then
    echo "notwritable"
  fi
  exit 0
fi

if [ -s ${prokkablastdb}/${refgenus} ] ; then
  datetime=$(date +%Y-%m-%d_%H-%M-%S)
  echo "Warning: found a previous Prokka reference database for the genus '${refgenus}; save it to '${prokkablastdb}/${refgenus}.backup_${datetime}'"
  mv -f ${prokkablastdb}/${refgenus} ${prokkablastdb}/${refgenus}.backup_${datetime}
fi
# add all the (representative) proteins in the dataset to the custom reference prot database for Prokka to search for similarities
ls ${inrefass}/*/*_genomic.gbff.gz > ${inrefass}_genomic_gbffgz_list
pgb2fdb=$(testpgb2fdb `head -n1 ${inrefass}_genomic_gbffgz_list`)
tmpd=${tmpdir}/mkprokkarefgendb
mkdir -p ${tmpd}
rm -f ${tmpd}/${refgenus}.faa.*  ${logdir}/prokka-genbank_to_fasta_db.log.*
if [ -s ${inrefass}_genomic_gbffgz_list ] ; then
  # extract protein sequences
  parallel --halt now,fail=1 -a ${inrefass}_genomic_gbffgz_list gzpgb2fdb "${pgb2fdb} {} ${tmpd}/${refgenus}.faa.{%} ${logdir}/prokka-genbank_to_fasta_db.log.{%}"
  if [ ${?} -gt 0 ] ; then
    >&2 echo "WARNING: could not build custom genus-specific BLAST db from files ${inrefass}/*/*_genomic.gbff ; Prokka will have to use default databases."
  else
    cat ${tmpd}/${refgenus}.faa.* > ${tmpd}/${refgenus}.faa && rm -f ${tmpd}/${refgenus}.faa.*
    ## cluster similar sequences
#    # use CD-HIT, 
#    # with thresholds: 90% seq id (-c 0.9), with seq identity defined as global (-G 1)
#	 # a mininum (unaligned!) length ratio of 80% (-s 0.8)
#	 # uses multiple threads (-T 0 <=> `nproc` threads), unlimited memory use (-M 0)
#    cdhit -i ${tmpd}/${refgenus}.faa -o ${tmpd}/${refgenus}_representative.faa -T 0 -M 0 -G 1 -s 0.8 -c 0.9 &> ${logdir}/make_prokka_ref_cdhit.log
#	 ndiscarded=$(grep -c 'Discarding invalid sequence' ${logdir}/cdhit.log)
#	 [[ ${ndiscarded} -gt 0 ]] && echo "Discarded ${ndiscarded} invalid sequence or sequence without identifier and description:" && grep -A2 'Discarding invalid sequence' | grep '>'
#    rm -fv ${tmpd}/${refgenus}.faa ${tmpd}/${refgenus}_representative.faa.clstr
#    cp -p ${tmpd}/${refgenus}_representative.faa ${prokkablastdb}/${refgenus}

    # use MMSeqs/Linclust,
    # with thresholds: 90% seq id (--min-seq-id 0.9), with seq identity defined over the (global) alignment length (--seq-id-mode 0)
    # a mininum (aligned!) length ratio of 80% (-c 0.8 --cov-mode 0)
	mmseqs easy-linclust ${tmpd}/${refgenus}.faa ${tmpd}/${refgenus}AnnotRef ${tmpd} &> ${logdir}/make_prokka_ref_mmseqs.log
    # replace database name for genus detection by Prokka 
    mv -f ${tmpd}/${refgenus}AnnotRef_rep_seq.fasta ${prokkablastdb}/${refgenus} && rm -rfv ${tmpd}/
    cd ${prokkablastdb}/
    makeblastdb -dbtype prot -in ${refgenus}
    cd -
  fi          
fi
