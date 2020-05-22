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

# parameters to be set:
#~ export recsamplesize=1000
if [ -z ${reccolid} ] ; then
 reccolid=1
fi
# derived parameters
if [ ${GeneRaxalgo} == 'reconciliation-samples' ] ; then
  export rectype='recsampling'
else
  export rectype='pointestimate'
fi
export reccol="generax_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_GeneRax_recs

gttorecdir=${coltreechains}/${collapsecond}/${replmethod}
grxlogs=${ptgdb}/logs/GeneRax
mkdir -p $grxlogs/${reccol}
export outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p ${outrecdir}

cd ${ptgtmp} 

#generaxcommonopt="-r UndatedDTL --max-spr-radius 5 --strategy SPR" # now a pipeline default

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  grlog=${grxlogs}/generax_global.log
  # use the same species tree file for every gene family, with no collapsed populations
  spetree=${speciestree}_clade_defs.nwk
  # this allows a single run of GeneRax, with built-in optimised load balance
  # generate a global family file i.e. job scheduling list and per-family parameter settings
  generaxfamfi=${alerec}/${reccol}_generax.families
  python ${ptgscripts}/make_generax_family_file.py --alignments ${cdsalifastacodedir} --out ${generaxfamfi}
  ${ptgscripts}/generax_global_mpi.sh ${generaxfamfi}

else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  spetree='Stree.nwk'
  # this dictate that every family need to be run independently, thus loosing the benefit of built-in optimised load balance
  
  # generate a family file i.e. parameter settings for each gene family
  generaxfamfidir=${alerec}/${reccol}_generax_families
  mkdir -p ${generaxfamfidir}/
  gttorecdir=${coltreechains}/${collapsecond}/${replmethod}
  python ${ptgscripts}/make_generax_family_file.py --perfam --alignments ${gttorecdir} --gene-trees ${gttorecdir} --out ${generaxfamfidir}
  
  tasklist=${generaxfamfidir}_list
  if [ -z ${genefamlist} ] ; then
    ${ptgscripts}/lsfullpath.py "${gttorecdir}/*_generax.families" > ${tasklist}
  else
    rm -f ${tasklist}
    for fam in $(cut -f1 ${genefamlist}) ; do
      ls ${generaxfamfidir}/${fam}*_generax.families 2> /dev/null
    done > ${tasklist} 
  fi
  
  grlog=${grxlogs}/generax_perfam.log
  export ncpus=1
  parallel -j ${ptgthreads} ${ptgscripts}/generax_perfam.sh :::: ${tasklist}
		
fi
