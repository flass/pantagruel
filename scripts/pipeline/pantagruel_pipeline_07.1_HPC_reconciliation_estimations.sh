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
#~ export ALEalgo='ALEml_undated'
#~ export recsamplesize=1000
#~ export ALEsourcenote='program compiled from source code from of https://github.com/ssolo/ALE/commits/63f0a3c964074a15f61fd45156ab9e10b5dd45ef'
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

tasklist=${coltreechains}/${collapsecond}/${colmethod}_Gtrees_list
ls ${coltreechains}/${collapsecond}/${colmethod}/*-Gtrees.nwk > $tasklist
alelogs=${ptgdb}/logs/ALE
mkdir -p $alelogs/${reccol}
outrecdir=${recs}/${collapsecond}/${colmethod}/${reccol}
mkdir -p $outrecdir

Njob=`wc -l $tasklist | cut -f1 -d' '`
qsubvars="tasklist=$tasklist, resultdir=$outrecdir, spetree=Stree.nwk, nrecs=${recsamplesize}, alealgo=${ALEalgo}"
qsub -J 1-$Njob -N ${reccol} -l select=1:ncpus=1:mem=20gb,walltime=24:00:00 -o $alelogs/${reccol} -j oe -v "$qsubvars" ${ptgscripts}/ALE_array_PBS.qsub

export reccoldate=$(date +%Y-%m-%d)
echo -e "${reccolid}\t${reccoldate}" > ${genetrees}/reccol
