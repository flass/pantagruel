#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

checkptgversion
checkfoldersafe ${database}

##############################################
## 03. Create and Populate SQLite database
##############################################

cd ${database}
### create and populate SQLite database

## Genome schema: initiate and populate
step1="step 1: initiating and populating database (genome-related tables)" && echo ${step1}
${ptgscripts}/pantagruel_sqlitedb_genome.sh "${database}" "${sqldbname}" "${genomeinfo}/assembly_metadata" "${genomeinfo}/assembly_info" "${protali}" "${protfamseqs}.tab" "${protorfanclust}" "${cdsorfanclust}" "${straininfo}" "${gblikeass}" "${gp2ass}"
checkexec "something went wrong while initiating the database"

# dump reference table for translation of genome assembly names into short identifier codes (uing UniProt "5-letter" code when available).
sqlite3 ${sqldb} "select assembly_id, code from assemblies;" | sed -e 's/|/\t/g' > ${database}/genome_codes.tab
checkexec "something went worng while generating 'genome_codes.tab' file"
sqlite3 ${sqldb} "select code, organism, strain from assemblies;" | sed -e 's/|/\t/g' > ${database}/organism_codes.tab
checkexec "something went worng while generating 'organism_codes.tab' file"
# and for CDS names into code of type "SPECIES_CDSTAG" (ALE requires such a naming)
# here split the lines with '|' and remove GenBank CDS prefix that is always 'lcl|'
sqlite3 ${sqldb} "select genbank_cds_id, cds_code from coding_sequences;" | sed -e 's/lcl|//g'  | sed -e 's/|/\t/g' > ${database}/cds_codes.tab
checkexec "something went worng while generating 'cds_codes.tab' file"
echo "succefully completed ${step1}"

# translates the header of the alignment files
step2="step 2: translation of the header of the alignment files from locus tags/protein ids to genome codes"
echo "${step2}"
for mol in prot cds ; do
  step2x="translation of the header of the $mol alignment files"
  alifastacodedir=${protali}/full_${mol}fam_alignments_species_code
  mkdir -p ${alifastacodedir}
  eval "export ${mol}alifastacodedir=${alifastacodedir}"
  # use multiprocessing python script
  ${ptgscripts}/lsfullpath.py "${protali}/full_${mol}fam_alignments/*.aln" > ${protali}/full_${mol}fam_alignment_list
  ${ptgscripts}/genbank2code_fastaseqnames.py ${protali}/full_${mol}fam_alignment_list ${database}/cds_codes.tab ${alifastacodedir} > ${ptgdb}/logs/genbank2code_fastaseqnames.${mol}.log
  checkexec "something went wrong during ${step2x}" "succefully completed ${step2x}"
done
if [[ "${compress}" == 'on' ]] ; then
  step2z="compressing redundant folders"
  echo "Pantagruel compress option (-z) ON: $step2z"
  cd ${protali}/
  for daln in full_protfam_alignments full_cdsfam_alignments ; do
    tar -czf ${daln}.tar.gz ${daln} && rm -r ${daln} || echo "Warning: could not succesfully compress '${protali}/${daln}/'; keep the full folder as is"
  done
fi
echo "succefully completed ${step2}"

## Phylogeny schema: initiate
step3="step 3: initiating database (phylogeny-related tables)" && echo ${step3}
sqlite3 ${sqldb} < ${ptgscripts}/pantagruel_sqlitedb_phylogeny_initiate.sql
checkexec "something went wrong during ${step3}" "succefully completed ${step3}"
