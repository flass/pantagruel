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

if [ -z "${ptgthreads}" ] ; then
  export ptgthreads=$(nproc)
fi

###############################################
## 07. Gene tree / Species tree reconciliations
###############################################

######################################################
## 07.1 Infer Gene tree / Species tree reconciliations
######################################################

### perform reconciliations with GeneRax, meaning the gene tree toppology is jointly optimized with the gene family evolution scenario

# parameters to be set: defaults:
#~ export ALEversion='v0.4'
#~ export ALEalgo='ALEml'
#~ export recsamplesize=1000
if [ -z ${reccolid} ] ; then
 reccolid=1
fi
# derived parameters
if [ ${GeneRaxalgo} == 'reconciliation-samples' ] ; then
  export rectype='recsampling'
  generaxopt="--reconciliation-samples ${recsamplesize}"
else
  export rectype='pointestimate'
fi
export reccol="generax_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_GenRax_recs

gttorecdir=${coltreechains}/${collapsecond}/${replmethod}
alelogs=${ptgdb}/logs/GeneRax
mkdir -p $alelogs/${reccol}
outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p ${outrecdir}

cd ${ptgtmp} 

generaxcommonopt="-r UndatedDTL --max-spr-radius 5 --strategy SPR"

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  spetree=${speciestree}_clade_defs.nwk
  # this allows a single run of GeneRax, with built-in optimised load balance
  # generate the family file i.e. job scheduling list
  generaxfamfi=${alerec}/${reccol}_generax.families
  python ${ptgscripts}/make_generax_family_file.py --alignments ${cdsalifastacodedir} --out ${generaxfamfi}
  mpiexec -np ${ptgthreads} generax ${generaxcommonopt} -s ${spetree} -f ${generaxfamfi} -p ${outrecdir} ${generaxopt}

else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  spetree='Stree.nwk'
  # this dictate that every family ned to be run independently, thus loosing the benefit of built-in optimised load balance
  echo "Error: NOT SUPPORTED YET; exit now" ; exit 1
  
  tasklist=${alerec}/${collapsecond}_${replmethod}_Gtree_list
  if [ -z ${genefamlist} ] ; then
    ${ptgscripts}/lsfullpath.py "${gttorecdir}/*-Gtree.nwk" > ${tasklist}
  else
    rm -f ${tasklist}
    for fam in $(cut -f1 ${genefamlist}) ; do
      ls ${gttorecdir}/${fam}*-Gtree.nwk 2> /dev/null
    done > ${tasklist} 
  fi
fi
