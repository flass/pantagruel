#!/usr/bin/env bash

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

### perform reconciliations with ecceTERA

# parameters to be set: defaults:
if [ -z "${reccolid}" ] ; then
 reccolid=1
fi
# derived parameters
if [ "${ecceTERAalgo}" == 'amalgamate' ] ; then
  # work on gene tree samples
  export rectype='ecceTERA_amalgamatedGchain'
  inputtrees="${coltreechains}/${collapsecond}/${replmethod}/*-Gtrees.nwk"
  inputtreetag="Gtrees"
else
  # work on single ML gene tree with branch supports
  inputtrees="${bayesgenetrees}/${collapsecond}/*.con.tre"
  inputtreetag="contre"
  if [ "${ecceTERAalgo:0:8}" == 'collapse' ] ; then
    # ecceTERA uses branch supports to 'collapse' the unresolved clades (collapsing appraoch similar in principle but different from Pantagruel's')
    export rectype="ecceTERA_collapsedGconstree_${ecceTERAalgo##*_}"
  else
    export rectype='ecceTERA'
  fi
fi
export reccol="${chaintype}_${rectype}_${reccolid}"


export tasklist=${alerec}/${collapsecond}_${replmethod}_${inputtreetag}_list
if [ -z "${genefamlist}" ] ; then
  ${ptgscripts}/lsfullpath.py "${inputtrees}" > ${tasklist}
else
  for fam in $(cut -f1 ${genefamlist}) ; do
     ls ${coltreechains}/${collapsecond}/${replmethod}/${fam}*-Gtrees.nwk
  done > ${tasklist}
fi
teralogs=${ptgdb}/logs/ecceTERA
mkdir -p ${teralogs}/${reccol}
export outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}
mkdir -p ${outrecdir}

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # use the same species tree file for every gene family, with no collapsed populations
  export spetree=${speciestree}.lsd.nwk ${recsamplesize} ${ALEalgo}
else
  # use a dedicated species tree file for each gene family, with population collapsed in accordance to the gene tree
  export spetree=Stree.nwk
fi


if [ "${resumetask}" == 'true' ] ; then
  # resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfgs in $(cat ${tasklist}) ; do
    bng=$(basename ${nfgs})
    bnalerec=${bng}*.reconciliationsFile_canonical_symmetric.txt
    nGdone=0
    if [ -e ${recs}/${collapsecond}/${replmethod}/${reccol}/${bnalerec} ] ; then
      nGdone=$(( ${nGdone} + 1))
    else
      echo ${nfgs}
    fi
  done > ${tasklist}_resumetasklist
  if [ ${nGdone} -gt 0 ] ; then
    echo "resume from previous computations: ${nGdone} G trees already reconciled"
    export tasklist=${tasklist}_resumetasklist
  else
    rm -f ${tasklist}_resumetasklist
  fi
fi

if [[ -z "${terabin}" ]] ; then
  export terabin="$(dirname $(readlink -f $(command -v ecceTERA)))"
fi

qsubvars="tasklist, outrecdir, spetree, terabin"
if [ ! -z "${fwdenv}" ] ; then
  qsubvars="${qsubvars}, ${fwdenv}"
fi
if [ ! -z "${alebin}" ] ; then
  qsubvars="${qsubvars}, alebin"
fi
if [ ! -z "${watchmem}" ] ; then
  qsubvars="${qsubvars}, watchmem"
fi

Njob=`wc -l ${tasklist} | cut -f1 -d' '`
[ ! -z ${topindex} ] &&  [ ${Njob} -gt ${topindex} ] && Njob=${topindex}
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))

for jobrange in ${jobranges[@]} ; do
 dlogs=${teralogs}/${reccol}/eccetera_${dtag}_${jobrange}
 mkdir -p ${dlogs}/

 case "$hpctype" in
  'PBS') 
    subcmd="qsub -J 1-$Njob -N ${reccol} -l select=1:ncpus=${ncpus}:mem=${mem}${parallelflags},walltime=${wth}:00:00 -o ${teralogs}/${reccol} -j oe -v \"$qsubvars\" ${ptgscripts}/ecceTERA_array_PBS.qsub"
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
    ${ptgscripts}/ecceTERA_array_LSF.bsub"
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

terasourcenote=""
pathterabin=$(readlink -f ${terabin})
terasourcenote="using ecceTERA software compiled from source; $(ecceTERA | grep version); binaries found at ${pathterabin}"
echo -e "${reccolid}\t${reccoldate}\t${terasourcenote}\t${reccol}" > ${alerec}/reccol
echo -e "\n# Reconciliation collection details:"
cat ${alerec}/reccol

