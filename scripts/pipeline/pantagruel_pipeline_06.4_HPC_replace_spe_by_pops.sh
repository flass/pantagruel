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
\t$0 [--ncpus (8)] [--mem (96)] [--wth (24)] [--hpctype ('PBS')] [--chunksize (100)] ptg_config_file\n
or for this help mesage:\n
\t$0 -h\n
\n
\tOption details:\n
\t --ncpus:\tnumber of parallel threads (default: 8).\n
\t --mem:\tmax memory allowance, in gigabytes (GB) (default: 32).\n
\t --wth:\tmax walltime, in hours (default: 24).\n
\t --hpctype:\t'PBS' for Torque scheduling system (default), or 'LSF' for IBM schedulling system.\n
\t --chunksize:\tnumber of families to be processed in one submitted job (default: 3000)\n
\t --fwdenv:\tnames of current environment variables you want to be passed down to the submitted jobs\n
\t\t\t if several, separate by ',' and enclose in quotes, e.g.: 'VAR1, VAR2' \n
\t --parallelflags:\tadditional flags to be specified for worker nodes' resource requests\n
\t\t\t for instance, specifying parallelflags=\"mpiprocs=1:ompthreads=8\" (with hpctype='PBS')\n
\t\t\t will force the use of OpenMP multi-threading instead of default MPI parallelism.\n
\t\t\t Note that parameter declaration syntax is not standard and differ depending on your HPC system;\n
\t\t\t also the relevant parameter flags may be different depending on the HPC system config.\n
\t\t\t Get in touch with your system admins to know the relevant flags and syntax.\n
\n
Note that the python scripts underlying this task require some non-standard packages.\n
It is advised to use Ancaonda to set up a stable environment to run these scripts\n
across the worker nodes of your HPC system. 
The job-submitted wrapper script will first attempt to load a conda environment\n
named 'ptgenv', which can be created (just once) using:\n
\t\`\`\`conda create -n ptgenv python=2.7 pip\n
\t\ \ \ conda activate ptgenv\n
\t\ \ \ pip install scipy numpy biopython bcbio-gff Cython igraph psycopg2\`\`\`\n
\n
If not found, the script will attempt to load Anaconda and attached python modules
with the \`module load anaconda3/personal\` command (not standard, likely not available on your cluster).
"

usage (){
  echo -e ${usagemsg}
}

ARGS=`getopt --options "h" --longoptions "ncpus:,mem:,wth:,hpctype:,chunksize:,parallelflags:,fwdenv:" --name "pantagruel_pipeline_06.4_HPC_replace_spe_by_pops.sh" -- "$@"`

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
envsourcescript="$1"
source ${envsourcescript}

if [ -z "${ncpus}" ] ; then
  export ncpus=8
fi
if [ -z "${mem}" ] ; then
  mem=96
fi
if [ -z "${wth}" ] ; then
  wth=24
fi
if [ -z "${hpctype}" ] ; then
  hpctype='PBS'
fi
if [ -z "${chunksize}" ] ; then
  chunksize=100
fi
if [ ! -z "${parallelflags}" ] ; then
  [ "${hpctype}" == 'PBS' ] && parallelflags=":${parallelflags}"
  [ "${hpctype}" == 'LSF' ] && parallelflags=" span[${parallelflags}]"
fi

if [ ! -z "${fwdenv}" ] ; then
  for varname in `echo "${fwdenv}" | sed -e "s/ *, */\n/g"` ; do
   export "${varname}"
  done
fi

################################################################################
## 06.4 Convert format of Bayesian gene trees and replace species by populations
################################################################################

## choice between pipline branches OPTION A2 ("${chaintype}" == 'fullgenetree') or OPTION B2 ("${chaintype}" == 'collapsed')
## is made within the submission script ${ptgscripts}/replSpePopinGs_array_PBS.qsub or ${ptgscripts}/replSpePopinGs_array_LSF.bsub
if [ -z ${replacecolid} ] ; then
 export replacecolid=1
fi
mkdir -p ${coltreechains}/${collapsecond}

mkdir -p ${ptgdb}/logs/replspebypop
export tasklist=${bayesgenetrees}_${collapsecond}_mbrun1t_list
ls ${bayesgenetrees}/${collapsecond}/${famprefix}*/*run1.t > $tasklist
repllogd=${ptgdb}/logs/replspebypop
export repllogs=${repllogd}/replace_species_by_pop_in_gene_trees
replacecoldate=$(date +%Y-%m-%d)
echo -e "${replacecolid}\t${replacecoldate}" > ${genetrees}/replacecol

export dtag="$(date +'%Y%m%d-%H%M%S')"
if [ "${resumetask}" == 'true' ] ; then
  rm -f ${tasklist}_resume*
  # resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfrun1t in $(cat $tasklist) ; do
   bnrun1t=$(basename $nfrun1t)
   bnGtre=${bnrun1t/.mb.run1.t/-Gtrees.nwk}
   bnStre=${bnrun1t/.mb.run1.t/-Stree.nwk}
   if [[ ! -e ${coltreechains}/${collapsecond}/${replmethod}/${bnGtre} || ! -e ${coltreechains}/${collapsecond}/${replmethod}/${bnStre} ]] ; then
     echo ${nfrun1t}
   fi
  done > ${tasklist}_resume
  export tasklist=${tasklist}_resume
fi
# PBS-submitted parallel job
# divide run in small chunks o be run in different jobs
Njob=`wc -l ${tasklist} | cut -f1 -d' '`
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
rm -f ${tasklist}_${dtag}_taskchunks
for jobrange in ${jobranges[@]} ; do
  replrun="${dtag}_${jobrange}"
  tail -n +$(echo $jobrange | cut -d'-' -f1) ${tasklist} | head -n ${chunksize} > ${tasklist}_${replrun}
  echo ${tasklist}_${replrun} >> ${tasklist}_${dtag}_taskchunks
done
Nchunk=`wc -l ${tasklist}_${dtag}_taskchunks | cut -f1 -d' '`
if [ ${Nchunk} -gt 2 ] ; then
  # run as an array job
  arrayspec=" -J 1-${Nchunk}"
else
  arrayspec=""
fi
case "${hpctype}" in 
  PBS)
    [ ${Nchunk} -gt 2 ] && arrayspec=" -J 1-${Nchunk}"
    subcmd="qsub${arrayspec} -N replSpePopinGs -l select=1:ncpus=${ncpus}:mem=${mem}gb${parallelflags},walltime=${wth}:00:00 -o ${repllogd} -j oe -V ${ptgscripts}/replSpePopinGs_array_PBS.qsub"
    ;;
  LSF)
    if [ ${Nchunk} -gt 2 ] ; then
      arrayspec="[1-${Nchunk}]"
      arrayjobtag='.%J'
    fi
    if [ ${wth} -le 12 ] ; then
      bqueue='normal'
    elif [ ${wth} -le 48 ] ; then
      bqueue='long'
    else
     bqueue='basement'
    fi
    memmb=$((${mem} * 1024)) 
    nflog="${repllogd}/replSpePopinGs.%I${arrayjobtag}.o"
    [ -z "${parallelflags}" ] && parallelflags="span[hosts=1]"
    subcmd="bsub -J replSpePopinGs${arrayspec} -R \"select[mem=${memmb}] rusage[mem=${memmb}] ${parallelflags}\" \
            -n${ncpus} -M${memmb} -q ${bqueue} \
            -o ${nflog} -e ${nflog} -env 'all' ${ptgscripts}/replSpePopinGs_array_LSF.bsub"
    ;;
esac
echo "$subcmd"
eval "$subcmd"

