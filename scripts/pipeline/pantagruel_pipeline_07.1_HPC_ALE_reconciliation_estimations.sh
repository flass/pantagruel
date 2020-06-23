#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018
# environment variables to be passed on to interface script:
export hpcscript=$(basename ${0})
hpcscriptdir=$(dirname ${0})
export defncpus=1
export defmem=96
export defwth=24
export defhpctype='PBS'
export defchunksize=1000
export withpython='false'

source ${hpcscriptdir}/pantagruel_pipeline_HPC_common_interface.sh "${@}"

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
if [ "${ALEalgo}" == 'ALEml_undated' ] ; then
  export rectype='undat'
  tag='u'
else
  export rectype='dated'
  tag=''
fi
export reccol="ale_${chaintype}_${rectype}_${reccolid}"
export recs=${alerec}/${chaintype}_ALE_recs

export tasklist=${alerec}/${collapsecond}_${replmethod}_Gtrees_list

${ptgscripts}/lsfullpath.py "${coltreechains}/${collapsecond}/${replmethod}/*-Gtrees.nwk" > $tasklist
alelogs=${ptgdb}/logs/ALE
mkdir -p ${alelogs}/${reccol}
export outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p $outrecdir
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  export spetree=${speciestree}.lsd.nwk
else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  export spetree='Stree.nwk'
fi

if [ "${resumetask}" == 'true' ] ; then
  rm -f ${tasklist}_resume
  # resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfgs in $(cat $tasklist) ; do
    bng=$(basename $nfgs)
    [ ${spetree} == 'Stree.nwk' ] && aleoutSpref=${bng/Gtrees/Stree} || aleoutSpref=$(basename ${spetree})
    bnalerec=${aleoutSpref}_${bng}.ale.${tag}ml_rec
    if [[ ! -e ${recs}/${collapsecond}/${replmethod}/${reccol}/${bnalerec} ]] ; then
     echo ${nfgs}
   fi
  done > ${tasklist}_resume
  export tasklist=${tasklist}_resume
fi

export dtag="$(date +'%Y%m%d-%H%M%S')"
export alebin="$alebin"
export watchmem="$watchmem"

qsubvars="tasklist, outrecdir, spetree, recsamplesize, ALEalgo, alebin, watchmem"
[ ! -z "${fwdenv}" ] && qsubvars="${qsubvars}, ${fwdenv}"
[ ! -z "${modulefile}" ] && qsubvars="${qsubvars}, modulefile"

Njob=`wc -l $tasklist | cut -f1 -d' '`
[ ! -z ${topindex} ] &&  [ ${Njob} -gt ${topindex} ] && Njob=${topindex}
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))

for jobrange in ${jobranges[@]} ; do
 dlogs=${alelogs}/${reccol}/ale_${dtag}_${jobrange}
 mkdir -p ${dlogs}/
 
 case "$hpctype" in
  'PBS') 
    subcmd="qsub -J 1-$Njob -N ${reccol} -l select=1:ncpus=${ncpus}:mem=${mem}${parallelflags},walltime=${wth}:00:00 -o $alelogs/${reccol} -j oe -v \"$qsubvars\" ${ptgscripts}/ale_array_PBS.qsub"
    ;;
  'LSF')
    if [ ${wth} -le 12 ] ; then
      bqueue='normal'
    elif [ ${wth} -le 48 ] ; then
      bqueue='long'
    else
	  bqueue='basement'
    fi
    arrayspec="[$jobrange]"
	[ ! -z ${maxatonce} ] && arrayspec="${arrayspec}%${maxatonce}"
    memmb=$((${mem} * 1024))
    nflog="${dlogs}/${reccol}.%J.%I.o"
	subcmd="bsub -J \"${reccol}${arrayspec}\" -q ${bqueue} \
	-R \"select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]\" -n${ncpus} -M${memmb} \
	-o ${nflog} -e ${nflog} -env \"${qsubvars}\" \
    ${ptgscripts}/ale_array_LSF.bsub"
	;;
  *)
    echo "Error: high-performance computer system '${hpctype}' is not supported; exit now"
    exit 1
	;;
 esac
 echo "${subcmd}"
 eval "${subcmd}"
done

export reccoldate=$(date +%Y-%m-%d)

if [[ -z "${alebin}" ]] ; then
  alebin=$(command -v ${ALEalgo})
fi

if [[ ! -z "$(echo ${alebin} | grep docker)" ]] ; then
  ALEsourcenote="using ALE Docker image $(docker image ls | grep alesuite | awk '{print $1,$3}')"
else
  pathalebin=$(readlink -f "${alebin}")
  alerepo=${pathalebin%%ALE/*}ALE/
  if [ -d ${alerepo} ] ; then
	alesrcvers=$(cd ${alerepo} && git log | head -n 1 | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
	alesrcorig=$(cd ${alerepo} && git remote -v | grep fetch | awk '{ print $2 }' 2> /dev/null && cd - > /dev/null)
  fi
  if [ ! -z "${alesrcvers}" ] ; then
    ALEsourcenote="using ALE software compiled from source; code origin: ${alesrcorig}; version ${alesrcvers}"
  else
    ALEheader="$(${ALEalgo} -h | head -n 1)"
    if [ -z "${ALEheader}" ] ; then
      ALEheader="using ALE software"
    fi
    ALEsourcenote="${ALEheader#*${ALEalgo} } binaries found at '${pathalebin}'"
  fi
fi

echo -e "${reccolid}\t${reccoldate}\t${ALEsourcenote}\t${reccol}" > ${alerec}/reccol
echo -e "\n# Reconciliation collection details:"
cat ${alerec}/reccol
