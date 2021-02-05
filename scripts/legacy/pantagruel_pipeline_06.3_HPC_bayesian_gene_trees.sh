#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019
# environment variables to be passed on to interface script:

echo "This is the 'usingGeneRax' branch of Pantagruel pipeline; using this module that runs 'MrBayes' bayesin gene tree inference method is not be of any use in the downstream pipeline (notably reconciliation with 'GeneRax') and may result into data that are incompatible with downstream modules."
echo "Alternatively, you can use the 'master' branch of the pipeline  for other reconciliation methods, ALE or ecceTERA: https://github.com/flass/pantagruel/tree/master."

export hpcscript=$(basename ${0})
hpcscriptdir=$(dirname ${0})
export defncpus=8
export defmem=16
export defwth=24
export defhpctype='PBS'
export defchunksize=1000
export withpython='false'

source ${hpcscriptdir}/pantagruel_pipeline_HPC_common_interface.sh "${@}"

############################
## 06.3 Bayesian gene trees
############################

## prpepare MrBayes run as a job array

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  export nexusaln4chains=${colalinexuscodedir}/${collapsecond}
else
  export nexusaln4chains=${colalinexuscodedir}/${collapsecond}/collapsed_alns
fi
export outputdir=${bayesgenetrees}/${collapsecond}
mkdir -p ${outputdir}

nruns=2
ngen=2000000
nchains=$(( ${ncpus} / ${nruns} ))

export tasklist=${nexusaln4chains}_ali_list
rm -f ${tasklist}
${ptgscripts}/lsfullpath.py "${nexusaln4chains}/*nex" > ${tasklist}

dtag=$(date +"%Y-%m-%d-%H-%M-%S")
append=''

# determine the set of numbered gene family prefixes to make separate folders
# and breakdown the load of files per folder
awk -F'/' '{print $NF}' ${tasklist} | grep -o "${famprefix}C[0-9]\{${ndiggrpfam}\}" | sort -u > ${nexusaln4chains}_ali_numprefixes
for pref in  `cat ${nexusaln4chains}_ali_numprefixes` ; do
  mkdir -p ${outputdir}/${pref}/
done

if [ "${resumetask}" == 'true' ] ; then
  #~ # following lines are for resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfaln in $(cat $tasklist) ; do
   bnaln=$(basename $nfaln)
   bncontre=${bnaln/nex/mb.con.tre}
   pref=$(echo ${bncontre} | grep -o "${famprefix}C[0-9]\{${ndiggrpfam}\}")
   if [ ! -e ${outputdir}/${pref}/${bncontre} ] ; then
    echo ${nfaln}
   fi
  done > ${tasklist}_resumetask_${dtag}
  export tasklist=${tasklist}_resumetask_${dtag}
  export mbmcmcopt='append=yes'
else
  export mbmcmcopt=''
fi
export mbmcmcpopt="Nruns=${nruns} Ngen=${ngen} Nchains=${nchains}"
qsubvars="tasklist, outputdir, mbmcmcopt, mbmcmcpopt, famprefix, ndiggrpfam"
if [ ! -z "${fwdenv}" ] ; then
  qsubvars="${qsubvars}, ${fwdenv}"
fi

Njob=`wc -l ${tasklist} | cut -f1 -d' '`
[ ! -z ${topindex} ] &&  [ ${Njob} -gt ${topindex} ] && Njob=${topindex}

jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
for jobrange in ${jobranges[@]} ; do
  dlogs=${ptglogs}/mrbayes/${chaintype}_mrbayes_trees_${collapsecond}_${dtag}_${jobrange}
  mkdir -p ${dlogs}/
 
  case "$hpctype" in
    'PBS') 
      subcmd="qsub -J ${jobrange} -N mb_panterodb -l select=1:ncpus=${ncpus}:mem=${mem}gb -l walltime=${wth}:00:00 \
	  -o ${dlogs} -v \"${qsubvars}\" ${ptgscripts}/mrbayes_array_PBS.qsub"
	  ;;
	'LSF')
	  if [ ${wth} -le 12 ] ; then
	    bqueue='normal'
	  elif [ ${wth} -le 48 ] ; then
	    bqueue='long'
	  else
	    bqueue='basement'
	  fi
	  memmb=$((${mem} * 1024))
	  nflog="${dlogs}/mrbayes.%J.%I.o"
	  subcmd="bsub -J \"mb_panterodb[$jobrange]\" -q ${bqueue} \
	  -R \"select[mem>${memmb}] rusage[mem=${memmb}] span[ptile=${ncpus}]\" -n${ncpus} -M${memmb} \
	  -o ${nflog} -e ${nflog} -env \"${qsubvars}\" \
	  ${ptgscripts}/mrbayes_array_LSF.bsub"
	  ;;
	*)
	  echo "Error: high-performance computer system '${hpctype}' is not supported; exit now"
	  exit 1
	  ;;
  esac
  echo "${subcmd}"
  eval "${subcmd}"
done
