#!/usr/bin/env bash
if [ -z "${1}" ] ; then
  echo "Usage: compress_ptgdb.sh ptg_env_file"
  exit 1
fi
source ${1}

#~ subfolders2="${protali}/collapsed_cdsfam_alignments_species_code ${protali}/full_cdsfam_fasta ${nrprotali} ${cdsalifastacodedir} ${protalifastacodedir} ${cdsalifastacodedir/_species_code/} ${protalifastacodedir/_species_code/}"

topfolders="${indata} ${seqdb} ${protali} ${database} ${funcannot} ${coregenome} ${genetrees} ${alerec} ${orthogenes} ${comparerecs}"

cd ${ptgdb}
for d in ${topfolders} ; do
  bnd=$(basename ${d})
  if [ -d ${bnd} ] ; then
    tar -czf ${bnd}.tar.gz && rm -r ${bnd}/
  fi
done

