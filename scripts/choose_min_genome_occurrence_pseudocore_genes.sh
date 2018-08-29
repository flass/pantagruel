#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

export raproot=$1
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source ${envsourcescript}

unset pseudocoremingenomes
mkdir -p ${coregenome}/pseudo-coregenome_sets/
# have glimpse of (almost-)universal unicopy gene family distribution and select those intended for core-genome tree given pseudocoremingenomes threshold
let "t = ($ngenomes * 3 / 4)" ; let "u = $t - ($t%20)" ; seq $u 20 $ngenomes | cat > ${ptgtmp}/mingenom ; echo "0" >> ${ptgtmp}/mingenom
Rscript --vanilla --silent ${ptgscripts}/select_pseudocore_genefams.r ${protali}/full_families_genome_counts-noORFans.mat ${database}/genome_codes.tab ${coregenome}/pseudo-coregenome_sets < ${ptgtmp}/mingenom

mv ${ptgdb}/environ_pantagruel_${ptgdbname}.sh ${ptgdb}/environ_pantagruel_${ptgdbname}.sh0 && \
 sed -e "s#'REPLACEpseudocoremingenomes'#$pseudocoremingenomes#" ${ptgdb}/environ_pantagruel_${ptgdbname}.sh0 > ${ptgdb}/environ_pantagruel_${ptgdbname}.sh && \
 rm ${ptgdb}/environ_pantagruel_${ptgdbname}.sh0
