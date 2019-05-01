#!/bin/bash

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

######################################################
## 07.1 Infer gene tree / Species tree reconciliations
######################################################

mkdir -p ${alerec}

### perform reconciliations with ALE

# parameters to be set: defaults:
#~ export ALEversion='v0.4'
#~ export ALEalgo='ALEml'
#~ export recsamplesize=1000
if [ -z ${reccolid} ] ; then
 reccolid=1
fi
# derived parameters
if [ ${ALEalgo} == 'ALEml_undated' ] ; then
  export rectype='undat'
else
  export rectype='dated'
fi
export reccol="ale_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_recs

tasklist=${alerec}/${collapsecond}_${replmethod}_Gtrees_list
ls ${coltreechains}/${collapsecond}/${replmethod}/*-Gtrees.nwk > $tasklist
alelogs=${ptgdb}/logs/ALE
mkdir -p $alelogs/${reccol}
outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p $outrecdir
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  spetree=${speciestree}.lsd.nwk ${recsamplesize} ${ALEalgo}
else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  spetree=Stree.nwk
fi

Njob=`wc -l $tasklist | cut -f1 -d' '`
qsubvars="tasklist=$tasklist, resultdir=$outrecdir, spetree=${spetree}, nrecs=${recsamplesize}, alealgo=${ALEalgo}"
qsub -J 1-$Njob -N ${reccol} -l select=1:ncpus=1:mem=20gb,walltime=24:00:00 -o $alelogs/${reccol} -j oe -v "$qsubvars" ${ptgscripts}/ale_array_PBS.qsub

export reccoldate=$(date +%Y-%m-%d)

if [ ! -z "$alebin" ] ; then
  alerepo=${alebin%%ALE/*}ALE/
  alesrcvers=$(cd ${alerepo} && git log | head -n 1 | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
  alesrcorig=$(cd ${alerepo} && git remote -v | grep fetch | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
  if [ ! -z "$alesrcvers" ] ; then
    ALEsourcenote="using ALE software compiled from source; code origin: ${alesrcorig}; version ${alesrcvers}"
  else
    ALEsourcenote="using ALE software binaries found at $(readlink -f ${alebin})"
else
  ALEsourcenote="using ALE Docker image $(docker image ls | grep alesuite | awk '{print $1,$3}')"
fi
echo -e "${reccolid}\t${reccoldate}\t${ALEsourcenote}" > ${alerec}/reccol
