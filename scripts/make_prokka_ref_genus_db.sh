#!/bin/bash


inrefass=${1}
refgenus=${2}
if [ -z ${3} ] ; then
  tmpdir=${3}
else
  tmpdir=$PWD
fiif [ -z ${4} ] ; then
  logdir=${4}
else
  logdir=$PWD
fi

prokkabin=$(which prokka)
prokkadir=$(dirname $(dirname $(readlink -f $prokkabin)))
prokkablastdb=${prokkadir}/db/genus

if [ -e ${prokkablastdb}/${refgenus} ] ; then
  datetime=$(date +%Y-%m-%d_%H-%M-%S)
  echo "Warning: found a previous Prokka reference database for the genus '${refgenus}; save it to '${prokkablastdb}/${refgenus}.backup_${datetime}'"
  mv -f ${prokkablastdb}/${refgenus} ${prokkablastdb}/${refgenus}.backup_${datetime}
fi
# add all the (representative) proteins in the dataset to the custom reference prot database for Prokka to search for similarities
ls ${inrefass}/*/*_genomic.gbff.gz > ${inrefass}_genomic_gbffgz_list
if [ -s ${inrefass}_genomic_gbffgz_list ] ; then
  parallel -a ${inrefass}_genomic_gbffgz_list 'gunzip -k'
  # extract protein sequences
  pgb2fdb="prokka-genbank_to_fasta_db ${inrefass}/*/*_genomic.gbff 1> ${tmpdir}/${refgenus}.faa 2> ${logdir}/prokka-genbank_to_fasta_db.log"
  eval "${pgb2fdb}"
  if [ ${?} -gt 0 ] ; then
    # try using the native perl
    /usr/bin/perl $prokkadir/bin/${pgb2fdb}
  fi
  if [ ${?} -gt 0 ] ; then
    >&2 echo "WARNING: could not build custom genus-specific BLAST db from files ${inrefass}/*/*_genomic.gbff ; Prokka will have to use default databases."
  else
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
