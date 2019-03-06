#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file [pseudocoremingenomes (int)]" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}
shift
if [ ! -z "${@}" ] ; then
  pseudocoremingenomes="${@}"
fi
mkdir -p ${coregenome}/pseudo-coregenome_sets/
# have glimpse of (almost-)universal unicopy gene family distribution and select those intended for core-genome tree given pseudocoremingenomes threshold
#~ let "t = ($ngenomes * 9 / 10)" ; let "u = $t - ($t%20)" ; seq $u 10 $ngenomes | sort -r > ${ptgtmp}/mingenom ; echo "0" >> ${ptgtmp}/mingenom

if [ ! -z "$pseudocoremingenomes" ] ; then
  if [[ "${pseudocoremingenomes}" == 'REPLACEpseudocoremingenomes' || "${pseudocoremingenomes}" == "${ngenomes}" ]] ; then
    echo "WARNING: Will rely on strict core genome definition to compute reference tree. This is often not advisable as the strict core genome can be very small." 1>&2
    echo "To choose a sensible 'pseudocore genomes' gene set, please run interactively '$ptgscripts/choose_min_genome_occurrence_pseudocore_genes.sh'." 1>&2
    export pseudocoremingenomes=${ngenomes}
  fi 
  # override interactivity with input file
  rm -f ${ptgtmp}/mingenom
  # if several values entered, iterate over them
  for p in ${pseudocoremingenomes} ; do
    echo ${p} >> ${ptgtmp}/mingenom
  done
  # repeat last value and thus make it the chosen value (other values will just have had their gene family set and heatmap representation computed for consultation)
  echo ${p} >> ${ptgtmp}/mingenom
  Rscript --vanilla --silent ${ptgscripts}/select_pseudocore_genefams.r \
   ${protali}/full_families_genome_counts-noORFans.mat ${database}/genome_codes.tab ${coregenome}/pseudo-coregenome_sets < ${ptgtmp}/mingenom
else
  # interactive call
  Rscript --vanilla --silent ${ptgscripts}/select_pseudocore_genefams.r \
   ${protali}/full_families_genome_counts-noORFans.mat ${database}/genome_codes.tab ${coregenome}/pseudo-coregenome_sets 2> ${ptgtmp}/set_pseudocoremingenomes
  eval "$(cat $ptgtmp/set_pseudocoremingenomes)"
fi
echo "set min number of genomes for inclusion in pseudo-core gene set as $pseudocoremingenomes"
mv ${envsourcescript} ${envsourcescript}0 && \
 sed -e "s#'REPLACEpseudocoremingenomes'#$pseudocoremingenomes#" ${envsourcescript}0 > ${envsourcescript} && \
 rm ${envsourcescript}0
echo "'export pseudocoremingenomes=$pseudocoremingenomes' recorded in ${envsourcescript}"
