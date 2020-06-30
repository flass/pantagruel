#!/bin/bash


#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (fl4@sanger.ac.uk), 13 January 2020
# environment variables to be passed on to interface script:
export hpcscript=$(basename ${0})
hpcscriptdir=$(dirname ${0})
export defncpus=1
export defmem=16
export defwth=24
export defhpctype='LSF'
export defchunksize=1000
export withpython='false'

source ${hpcscriptdir}/pantagruel_pipeline_HPC_common_interface.sh "${@}"

######################################################
## 07.1 Infer gene tree / Species tree reconciliations
######################################################

### perform reconciliations with GeneRax, meaning the gene tree toppology is jointly optimized with the gene family evolution scenario

# parameters to be set:
#~ export recsamplesize=1000
if [ -z "${reccolid}" ] ; then
 reccolid=1
fi
# derived parameters
if [ "${GeneRaxalgo}" == 'reconciliation-samples' ] ; then
  export rectype='recsampling'
else
  export rectype='pointestimate'
fi
export reccol="generax_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_GeneRax_recs

gttorecdir=${coltreechains}/${collapsecond}/${replmethod}
grxlogs=${ptgdb}/logs/GeneRax
mkdir -p ${grxlogs}/${reccol}
export outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p ${outrecdir}

cd ${ptgtmp} 

# load GeneRax support modules
if [ ! -z "${modulefile}" ] ; then
  source ${modulefile}
fi


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
#generaxcommonopt="-r UndatedDTL --max-spr-radius 5 --strategy SPR" # now a pipeline default

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  export spetree=${speciestree}_clade_defs.nwk
  export alitorecdir="${cdsalifastacodedir} --skip-abs-gt"
else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  export spetree='Stree.nwk'
  export spetreedir=${gttorecdir}
  export alitorecdir=${gttorecdir}
  # this dictate that every family need to be run independently, thus loosing the benefit of built-in optimised load balance
fi
if [[ "${chaintype}" == 'fullgenetree' && "${GeneRaxalgo}" =~ 'global' ]] ; then
  # using the same species tree allows a single run of GeneRax, with built-in optimised load balance
  echo "detected 'global' keyword in reconciliation algorithm"
  echo "will run GeneRax on the whole pangenome with global parameter estimation"
  # generate a global family file i.e. job scheduling list and per-family parameter settings
  generaxfamfi=${alerec}/${reccol}_generax.families
  step1="create a family file i.e. parameter settings for the whole pangenome gene family set"
  python2.7 ${ptgscripts}/make_generax_family_file.py --alignments ${alitorecdir} --out ${generaxfamfi} --gene-trees ${gttorecdir}
  checkexec "failed to ${step1}" "successfully ${step1/create/created}"
  ls ${generaxfamfi}
  step2="run GeneRax on all pangenome genes at once"
  
  echo "submitting MPI job to ${step2} (using GeneRax built-in optimised load balance on ${ncpus} cores)"
  case "$hpctype" in
    'PBS') 
        subcmd="qsub -l walltime=${wth}:00:00 -l select=1:ncpus=${ncpus}:mem=${mem}gb -N \"${reccol}\" \
	    -o ${grxlogs}/${reccol}/generax_global_mpi.log -j oe -V ${ptgscripts}/generax_global_mpi.qsub ${generaxfamfi}"
        ;;
    'LSF')
	  bqueue='parallel'
      memmb=$((${mem} * 1024))
	  nflog="${grxlogs}/${reccol}/generax_global_mpi.%J.log"
      subcmd="bsub -q parallel -J \"${reccol}\" -R \"select[mem>${memmb}] rusage[mem=${memmb}] span[ptile=${ncpuperhost}]\" -M ${memmb} -n ${ncpus} \
       -o ${grlog} -e ${grlog} -env 'all' ${ptgscripts}/generax_global_mpi.bsub ${generaxfamfi}"
	  ;;
    *)
	  echo "Error: high-performance computer system '${hpctype}' is not supported; exit now"
      exit 1
	  ;;
  esac
  echo "${subcmd}"
  eval "${subcmd}"

else
  echo "will run GeneRax independently on each pangenome gene family, with family-specific parameter estimation"
  # generate a family file i.e. parameter settings for each gene family
  step1="create a family file i.e. parameter settings for each gene family"
  generaxfamfidir=${alerec}/${reccol}_generax_families
  mkdir -p ${generaxfamfidir}/
  # collapsed-replaced alignments are in the same folder as collapsed-replaced gene trees
  python2.7 ${ptgscripts}/make_generax_family_file.py --per-family --alignments ${alitorecdir} \
   --gene-trees ${gttorecdir} --out ${generaxfamfidir} --gftag '.generax_families'
  checkexec "failed to ${step1}" "successfully ${step1/create/created}"
  ls -d ${generaxfamfidir}
  
  tasklist=${generaxfamfidir}_list
  if [ -z "${genefamlist}" ] ; then
    ${ptgscripts}/lsfullpath.py "${generaxfamfidir}/*.generax_families" > ${tasklist}
  else
    rm -f ${tasklist}
    for fam in $(cut -f1 ${genefamlist}) ; do
      ls ${generaxfamfidir}/${fam}*.generax_families 2> /dev/null
    done > ${tasklist} 
  fi
  
  Njob=`wc -l ${tasklist} | cut -f1 -d' '`
  [ ! -z "${topindex}" ] &&  [ ${Njob} -gt ${topindex} ] && Njob=${topindex}
  jobranges=($(${ptgscripts}/get_jobranges.py ${chunksize} ${Njob}))

  step2="run GeneRax on each pangenome gene family in parallel"
  echo "submitting array job(s) to ${step2} (using ${ncpus} thread(s) per process) with the following range: ${jobranges[@]}"
  for jobrange in ${jobranges[@]} ; do
    dlogs=${grxlogs}/${reccol}/generax_perfam_array_${dtag}_${jobrange}
    mkdir -p ${dlogs}/
    case "${hpctype}" in
      'PBS') 
        subcmd="qsub -J ${jobrange} -l walltime=${wth}:00:00 -l select=1:ncpus=${ncpus}:mem=${mem}gb -N \"${reccol}\" \
	    -o ${dlogs} -j oe -V ${ptgscripts}/generax_perfam_array.qsub"
        ;;
      'LSF')
	    bqueue='parallel'
        arrayspec="[${jobrange}]"
	    [ ! -z "${maxatonce}" ] && arrayspec="${arrayspec}%${maxatonce}"
        memmb=$((${mem} * 1024))
        nflog="${dlogs}/generax_perfam.%J.%I.log"
  
        subcmd="bsub -J \"${reccol}${arrayspec}\" -q ${bqueue} -R \"select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]\" \
		 -M ${memmb} -n ${ncpus} -o ${nflog} -e ${nflog} -env 'all' ${ptgscripts}/generax_perfam_array.bsub"
	    ;;
	  *)
	    echo "Error: high-performance computer system '${hpctype}' is not supported; exit now"
        exit 1
	    ;;
    esac
    echo "${subcmd}"
    eval "${subcmd}"
  done
  
fi

echo -e "${reccolid}\t${reccoldate}\t${GRsourcenote}\t${reccol}" > ${alerec}/reccol
echo -e "\n# Reconciliation collection details:"
cat ${alerec}/reccol
