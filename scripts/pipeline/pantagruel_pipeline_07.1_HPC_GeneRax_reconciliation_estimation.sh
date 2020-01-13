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
export defmem=64
export defwth=168
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

  case "$hpctype" in
    'LSF')
	  bqueue='parallel'
      memmb=$((${mem} * 1024))
	  
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
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  spetree='Stree.nwk'
  # this dictate that every family ned to be run independently, thus loosing the benefit of built-in optimised load balance
  
  # generate a family file i.e. parameter settings for each gene family
  generaxfamfidir=${alerec}/${reccol}_generax_families
  mkdir -p ${generaxfamfidir}/
  gttorecdir=${coltreechains}/${collapsecond}/${replmethod}
  python ${ptgscripts}/make_generax_family_file.py --perfam --alignments ${gttorecdir} --gene-trees ${gttorecdir} --out ${generaxfamfidir}
  
  export tasklist=${generaxfamfidir}_list
  if [ -z ${genefamlist} ] ; then
    ${ptgscripts}/lsfullpath.py "${generaxfamfidir}/*.generax_families" > ${tasklist}
  else
    rm -f ${tasklist}
    for fam in $(cut -f1 ${genefamlist}) ; do
      ls ${generaxfamfidir}/${fam}*.generax_families 2> /dev/null
    done > ${tasklist} 
  fi
  
  grlog=${grxlogs}/generax_perfam.log
  
  Njob=`wc -l ${tasklist} | cut -f1 -d' '`
  [ ! -z ${topindex} ] &&  [ ${Njob} -gt ${topindex} ] && Njob=${topindex}
  jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))

  for jobrange in ${jobranges[@]} ; do
    dlogs=${grxlogs}/${reccol}/generax_perfam_${dtag}_${jobrange}.log
    mkdir -p ${dlogs}/
    case "$hpctype" in
      'LSF')
	    bqueue='parallel'
        arrayspec="[$jobrange]"
	    [ ! -z ${maxatonce} ] && arrayspec="${arrayspec}%${maxatonce}"
        memmb=$((${mem} * 1024))
        nflog="${dlogs}/${reccol}.%J.%I.o"
  
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
