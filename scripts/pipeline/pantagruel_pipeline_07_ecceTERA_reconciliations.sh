#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file [gene_fam_list]" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

if [ ! -z "$2" ] ; then
  export genefamlist="$2"
fi

checkptgversion
checkfoldersafe ${alerec}

###############################################
## 07. Gene tree / Species tree reconciliations    with ecceTERA
###############################################

######################################################
## 07.1 Infer gene tree / Species tree reconciliations
######################################################

### perform reconciliations with ecceTERA

ecceTERA species.file=BRADYC021242-collapsed-Stree.nwk gene.file=BRADYC021242-collapsed-Gtrees.nwk verbose=true amalgamate=true print.newick=true print.info=true &> BRADYC021242-collapsed-Gtrees.nwk.ecceTERA_withcollapsedStree_amalgamate.log
mv geneTree BRADYC021242-collapsed-Gtrees.nwk.ecceTERA_withcollapsedStree_amalgamate.geneTree
mv speciesTree BRADYC021242-collapsed-Gtrees.nwk.ecceTERA_withcollapsedStree_amalgamate.speciesTree

# parameters to be set: defaults:
#~ export ALEversion='v0.4'
#~ export ALEalgo='ALEml'
#~ export recsamplesize=1000
if [ -z ${reccolid} ] ; then
 reccolid=1
fi
# derived parameters
if [ ${chaintype} == 'collapsed' ] ; then
  # work on gene tree samples
  export rectype='ecceTERA_amalgamate'
else
  # work on single ML gene tree with branch supports, used to 'collapse' (in the [much similar] ecceTERA meaning) the unresolved clades
  export rectype='ecceTERA_collapse1'
fi
export reccol="ale_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_ecceTERA_recs

tasklist=${alerec}/${collapsecond}_${replmethod}_Gtrees_list
if [ -z ${genefamlist} ] ; then
  ${ptgscripts}/lsfullpath.py "${coltreechains}/${collapsecond}/${replmethod}/*-Gtrees.nwk" > ${tasklist}
else
  rm -f ${tasklist}
  for fam in $(cut -f1 ${genefamlist}) ; do
    ls ${coltreechains}/${collapsecond}/${replmethod}/${fam}*-Gtrees.nwk >> ${tasklist}
  done
fi
alelogs=${ptgdb}/logs/ALE
mkdir -p $alelogs/${reccol}
outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p $outrecdir

cd ${ptgtmp} 
