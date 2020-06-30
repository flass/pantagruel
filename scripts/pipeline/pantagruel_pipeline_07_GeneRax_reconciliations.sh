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
if [[ "${GeneRaxalgo}" =~ 'reconciliation-samples' ]] ; then
  export rectype='recsampling'
else
  export rectype='pointestimate'
fi
export reccol="generax_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_GeneRax_recs
export reccoldate=$(date +%Y-%m-%d)

gttorecdir=${coltreechains}/${collapsecond}/${replmethod}
grxlogs=${ptgdb}/logs/GeneRax
mkdir -p ${grxlogs}/${reccol}
export outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p ${outrecdir}

# recording the software version that was used
if [[ -z "${grbin}" ]] ; then
  grbin=$(command -v generax)
fi
if [[ -z "${grbin}" ]] ; then
  echo "Error: could not find command 'generax' in the path; please provide the path to the generax executable file through env variable \$grbin; exit now"
  exit 1
else
  GRheader="$(${grbin} -h | grep '\[00:00:' | awk '{ print $2,$3 }')"
  if [ -z "${GRheader}" ] ; then
    GRheader="GeneRax"
  fi
  GRvtag=$(echo ${GRheader} | awk '{ print $NF }')
  pathgrbin=$(readlink -f "${grbin}")
  grrepo=${pathgrbin%%GeneRax/*}GeneRax/
  if [ -d ${grrepo} ] ; then
	grsrcvers=$(cd ${grrepo} && git log | head -n 1 | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
	grsrcorig=$(cd ${grrepo} && git remote -v | grep fetch | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
  fi
  if [ ! -z "${grsrcvers}" ] ; then
    GRsourcenote="using GeneRax software (version ${GRvtag}) compiled from source; code origin: ${grsrcorig}; code version ${grsrcvers}"
  else
    GRsourcenote="${GRheader} binaries found at '${pathgrbin}'"
  fi
fi

# generaxcommonopt="-r UndatedDTL --max-spr-radius 5 --strategy SPR" # now a pipeline default

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  export spetree=${speciestree}_clade_defs.nwk
  # expect that task 6 has been skipped altogether; proceed with the gene family selection
  echo "generate list of gene familes to be reconciled"
  allfamlist=${alerec}/cdsfams_minsize4
  python2.7 ${ptgscripts}/pantagruel_sqlitedb_query_gene_fam_sets.py --db=${sqldb} --outprefix='cdsfams' --dirout=${alerec} --base.query="${basequery}" \
   --famsets.min.sizes="4"
  checkexec "could not generate gene family list ; exit now" "succesfully generated gene family list : $allfamlist"
  if [ -z "${genefamlist}" ] ; then
    famlist=${allfamlist}
  else
    echo "will restrict the gene families to be processed to those listed in '${genefamlist}':"
    cut -f1 ${genefamlist} | head
    ngf=$(wc -l ${genefamlist} | cut -d' ' -f1)
    if [ $ngf -gt 6 ] ; then
      echo "... ($ngf families)"
    fi
    bngenefamlist=$(basename ${genefamlist})
    famlist=${genetrees}/${bngenefamlist}_cdsfams_minsize4
    for fam in $(cut -f1 ${genefamlist}) ; do
      grep ${fam} ${allfamlist}
    done > ${famlist}
    nrestrfam=$(wc -l $famlist)
    checkexec "could not generate restricted gene family list ; exit now" "succesfully generated restricted gene family lists :\n$famlist (containing $(echo $nrestrfam | cut -d' ' -f1) gene families)"
  fi
  alntasklist=${famlist}_aln_list
  rm -f ${alntasklist}*
  for fam in $(cut -f1 ${famlist}) ; do
    ls ${cdsalifastacodedir}/${fam}.codes.aln >> ${alntasklist}
  done
  export mkgrfamfiopts="--alignment_list ${alntasklist}"
else
  # use a pre-computed gene tree with collapsed rake clades and a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  export spetree='Stree.nwk'
  export spetreedir=${gttorecdir}
  # this dictate that every family need to be run independently, thus loosing the benefit of built-in optimised load balance
  export mkgrfamfiopts="--alignments ${gttorecdir} --gene-trees ${gttorecdir}"
fi
if [[ "${chaintype}" == 'fullgenetree' && "${GeneRaxalgo}" =~ 'global' ]] ; then
  # using the same species tree allows a single run of GeneRax, with built-in optimised load balance
  echo "detected 'global' keyword in reconciliation algorithm"
  echo "will run GeneRax on the whole pangenome with global parameter estimation"
  # generate a global family file i.e. job scheduling list and per-family parameter settings
  generaxfamfi=${alerec}/${reccol}_generax.families
  step1="create a family file i.e. parameter settings for the whole pangenome gene family set"
  python2.7 ${ptgscripts}/make_generax_family_file.py ${mkgrfamfiopts} --out ${generaxfamfi}
  checkexec "failed to ${step1}" "successfully ${step1/create/created}"
  step2="run GeneRax on all pangenome genes at once"
#  grlog=${grxlogs}/generax_global.log
  export ncpus=${ptgthreads}
  echo "${step2} (using GeneRax built-in optimised load balance on ${ncpus} cores)"
  ${ptgscripts}/generax_global_mpi.sh ${generaxfamfi}
  checkexec "failed to ${step2}" "successfully ${step2/run/ran}"
else  
  echo "will run GeneRax independently on each pangenome gene family, with family-specific parameter estimation"
  # generate a family file i.e. parameter settings for each gene family
  step1="create a family file i.e. parameter settings for each gene family"
  generaxfamfidir=${alerec}/${reccol}_generax_families
  mkdir -p ${generaxfamfidir}/
  # collapsed-replaced alignments are in the same folder as collapsed-replaced gene trees
  python2.7 ${ptgscripts}/make_generax_family_file.py ${mkgrfamfiopts} --per-family --out ${generaxfamfidir} --gftag '.generax_families'
  checkexec "failed to ${step1}" "successfully ${step1/create/created}"
  
  tasklist=${generaxfamfidir}_list
  if [ -z "${genefamlist}" ] ; then
    ${ptgscripts}/lsfullpath.py "${generaxfamfidir}/*.generax_families" > ${tasklist}
  else
    rm -f ${tasklist}
    for fam in $(cut -f1 ${genefamlist}) ; do
      ls ${generaxfamfidir}/${fam}*.generax_families 2> /dev/null
    done > ${tasklist} 
  fi
  
#  grlog=${grxlogs}/generax_perfam.log
  export ncpus=1
  step2="run GeneRax on each pangenome gene family in parallel"
  echo "${step2} (using GNU parallel on ${ptgthreads} cores, ${ncpus} thread(s) per process)"
  parallel -j ${ptgthreads} --results ${grxlogs}/${reccol}/generax_perfam.{/.}.log ${ptgscripts}/generax_perfam.sh :::: ${tasklist}
  checkexec "failed to ${step2}" "successfully ${step2/run/ran}"
		
fi

echo -e "${reccolid}\t${reccoldate}\t${GRsourcenote}\t${reccol}" > ${alerec}/reccol
echo -e "\n# Reconciliation collection details:"
cat ${alerec}/reccol
