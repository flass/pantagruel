#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

export raproot=$1
envsourcescript=${raproot}/environ_pantagruel.sh
source $envsourcescript

unset pseudocoremingenomes
mkdir -p ${coregenome}/pseudo-coregenome_sets/
# have glimpse of (almost-)universal unicopy gene family distribution and select those intended for core-genome tree given pseudocoremingenomes threshold
let "t = ($ngenomes * 3 / 4)" ; let "u = $t - ($t%20)" ; seq $u 20 $ngenomes | cat > ${raptmp}/mingenom ; echo "0" >> ${raptmp}/mingenom
Rscript --vanilla --silent ${ptgscripts}/select_pseudocore_genefams.r ${protali}/full_families_genome_counts-noORFans.mat ${database}/genome_codes.tab ${coregenome}/pseudo-coregenome_sets < ${raptmp}/mingenom

mv ${raproot}/environ_pantagruel.sh ${raproot}/environ_pantagruel.sh0 && \
 sed -e "s#'REPLACEpseudocoremingenomes'#$pseudocoremingenomes#" ${raproot}/environ_pantagruel.sh0 > ${raproot}/environ_pantagruel.sh
 
