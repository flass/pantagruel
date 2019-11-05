#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019

usagemsg="
Usage:\n
 For running the script: with [option (default)]\n
\t$0 [--ncpus (8)] [--mem (32)] [--wth (24)] [--hpctype ('PBS')] [--chunksize (3000)] ptg_config_file\n
or for this help mesage:\n
\t$0 -h\n
\n
\tOption details:\n
\t --ncpus:\tnumber of parallel threads (default: 8).\n
\t --mem:\tmax memory allowance, in gigabytes (GB) (default: 32).\n
\t --wth:\tmax walltime, in hours (default: 24).\n
\t --hpctype:\t'PBS' for Torque scheduling system (default), or 'LSF' for IBM schedulling system.\n
\t --chunksize:\tnumber of families to be processed in one submitted job (default: 1000)\n
\t --fwdenv:\tnames of current environment variables you want to be passed down to the submitted jobs\n
\t\t\t if several, separate by ',' and enclose in quotes, e.g.: 'VAR1, VAR2' \n
\t --parallelflags:\tadditional flags to be specified for worker nodes' resource requests\n
\t\t\t for instance, specifying parallelflags=\"mpiprocs=1:ompthreads=8\" (with hpctype='PBS')\n
\t\t\t will force the use of OpenMP multi-threading instead of default MPI parallelism.\n
\t\t\t Note that parameter declaration syntax is not standard and differ depending on your HPC system;\n
\t\t\t also the relevant parameter flags may be different depending on the HPC system config.\n
\t\t\t Get in touch with your system admins to know the relevant flags and syntax.\n"

usage (){
  echo -e ${usagemsg}
}

ARGS=`getopt --options "h" --longoptions "ncpus:,mem:,wth:,hpctype:,chunksize:,parallelflags:,fwdenv:" --name "pantagruel_pipeline_06.3_HPC_bayesian_gene_trees.sh" -- "$@"`

#Bad arguments
if [ ${?} -ne 0 ];
then
  usage ; exit 1
fi

eval set -- "${ARGS}"

while true;
do
  case "${1}" in
    -h|--help) 
      echo usage
      exit 0;;
	
    --ncpus)
      testmandatoryarg "${1}" "${2}"
      export ncpus="${2}"
      echo "will run jobs using ${ncpus} parallel threads"
      shift 2;;
	  
	--mem)
      testmandatoryarg "${1}" "${2}"
      export mem="${2}"
      echo "will request ${mem} GB for submitted jobs"
      shift 2;;
	  
	--wth)
      testmandatoryarg "${1}" "${2}"
      export wth="${2}"
      echo "will request ${wth} hours walltime for submitted jobs"
      shift 2;;
	
	--chunksize)
      testmandatoryarg "${1}" "${2}"
      export chunksize="${2}"
      echo "will divide the work load into chunks of ${chunksize} individual tasks (gene trees)"
      shift 2;;
	
	--parallelflags)
      testmandatoryarg "${1}" "${2}"
      export parallelflags="${2}"
      echo "will add the following flags to the job submission resource request: '${parallelflags}'"
      shift 2;;
	
	--fwdenv)
      testmandatoryarg "${1}" "${2}"
      export fwdenv="${2}"
      echo "will forward the following environment variables to the submitted jobs' envronment: '${fwdenv}'"
      shift 2;;
	  
	--)
      shift
      break;;
    
    *)
	  echo "Error: this argument is not supported: ${1}"
	  exit 1;;
	  
  esac
done

if [ -z "${1}" ] ; then
 echo -e "Error: missing mandatory parameter: pantagruel config file\n"
 usage
 exit 1
fi
envsourcescript="${1}"
source ${envsourcescript}

if [ -z "${ncpus}" ] ; then
  export ncpus=8
fi
if [ -z "${mem}" ] ; then
  mem=24
fi
if [ -z "${wth}" ] ; then
  wth=24
fi
if [ -z "${hpctype}" ] ; then
  hpctype='PBS'
fi
if [ -z "${chunksize}" ] ; then
  chunksize=1000
fi
if [ ! -z "${parallelflags}" ] ; then
  [ "${hpctype}" == 'PBS' ] && parallelflags=":${parallelflags}"
  [ "${hpctype}" == 'LSF' ] && parallelflags=" span[${parallelflags}]"
fi

############################
## 06.3 Bayesian gene trees
############################

## run mrbayes on collapsed alignments
export nexusaln4chains=${colalinexuscodedir}/${collapsecond}/collapsed_alns
export mboutputdir=${bayesgenetrees}/${collapsecond}
mkdir -p ${mboutputdir}
nruns=2
nchains=$(( ${ncpus} / ${nruns} ))
tasklist=${nexusaln4chains}_ali_list
rm -f ${tasklist}
${ptgscripts}/lsfullpath.py "${nexusaln4chains}/*nex" > ${tasklist}
dtag=$(date +"%Y-%m-%d-%H-%M-%S")
append=''

# determine the set of numbered gene family prefixes to make separate folders
# and breakdown the load of files per folder
awk -F'/' '{print $NF}' ${mbtasklist} | grep -o "${famprefix}C[0-9]\{${ndiggrpfam}\}" | sort -u > ${nexusaln4chains}_ali_numprefixes
for pref in  `cat ${nexusaln4chains}_ali_numprefixes` ; do
  mkdir -p ${mboutputdir}/${pref}/
done

if [ "${resumetask}" == 'true' ] ; then
  #~ # following lines are for resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfaln in $(cat $tasklist) ; do
   bnaln=$(basename $nfaln)
   bncontre=${bnaln/nex/mb.con.tre}
   pref=$(echo ${bncontre} | grep -o "${famprefix}C[0-9]\{${ndiggrpfam}\}")
   if [ ! -e ${mboutputdir}/${pref}/${bncontre} ] ; then
    echo ${nfaln}
   fi
  done > ${tasklist}_resumetask_${dtag}
  tasklist=${tasklist}_resumetask_${dtag}
  export mbmcmcopt='append=yes'
else
  export mbmcmcopt=''
fi
export mbmcmcpopt="Nruns=${nruns} Ngen=2000000 Nchains=${nchains}"
qsubvars="tasklist=${tasklist}, outputdir=${mboutputdir}"
if [ ! -z "${fwdenv}" ] ; then
  qsubvars="${qsubvars}, ${fwdenv}"
fi

Njob=`wc -l ${tasklist} | cut -f1 -d' '`
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
for jobrange in ${jobranges[@]} ; do
 dlogs=${ptglogs}/mrbayes/${chaintype}_mrbayes_trees_${collapsecond}_${dtag}_${jobrange}
 mkdir -p ${dlogs}/
 
 case "$hpctype" in
    'PBS') 
      subcmd="qsub -J ${jobrange} -N mb_panterodb -l select=1:ncpus=${ncpus}:mem=${mem}gb -l walltime=${wth}:00:00 \
	  -o ${dlogs} -v \"$qsubvars, mbmcmcopt, mbmcmcpopt, famprefix\" ${ptgscripts}/mrbayes_array_PBS.qsub"
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
	  -o ${nflog} -e ${nflog} -env \"$qsubvars, mbmcmcopt, mbmcmcpopt, famprefix\" \
	  ${ptgscripts}/mrbayes_array_LSF.bsub"
	  ;;
	*)
	  echo "Error: high-performance computer system '$hpctype' is not supported; exit now"
	  exit 1;;
  esac
  echo "${subcmd}"
  eval "${subcmd}"
done
