#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 6 Nov 2019

# environment variables inherited from parent script:
#hpcscript
#defncpus
#defmem
#defwth
#defhpctype
#defchunksize
#withpython

usagemsg="
Usage:\n
 For running the script: with [option (default)]\n
\t${hpcscript} [--ncpus (${defncpus})] [--mem (${defmem})] [--wth (${defwth})] [--hpctype ('${defhpctype}')] [--chunksize (${defchunksize})] ptg_config_file\n
or for this help mesage:\n
\t${hpcscript} -h\n
\n
\tOption details:\n
\t --ncpus:\tnumber of parallel threads (default: ${defncpus}).\n
\t --mem:\tmax memory allowance, in gigabytes (GB) (default: ${defmem}).\n
\t --wth:\tmax walltime, in hours (default: ${defwth}).\n
\t --hpctype:\t'PBS' for Torque scheduling system, or 'LSF' for IBM schedulling system (default: '${defhpctype}').\n
\t --chunksize:\tnumber of families to be processed in one submitted job (default: ${defchunksize})\n
\t --fwdenv:\tnames of current environment variables you want to be passed down to the submitted jobs\n
\t\t\t if several, separate by ',' and enclose in quotes, e.g.: 'VAR1, VAR2' \n
\t --parallelflags:\tadditional flags to be specified for worker nodes' resource requests\n
\t\t\t for instance, specifying parallelflags=\"mpiprocs=1:ompthreads=8\" (with hpctype='PBS')\n
\t\t\t will force the use of OpenMP multi-threading instead of default MPI parallelism.\n
\t\t\t Note that parameter declaration syntax is not standard and differ depending on your HPC system;\n
\t\t\t also the relevant parameter flags may be different depending on the HPC system config.\n
\t\t\t Get in touch with your system admins to know the relevant flags and syntax.\n
"
pythonmsg="
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
  [ "${withpython}" == 'true' ] && echo -e ${pythonmsg}
}

# function to check the presence of value in parsing option arguments
testmandatoryarg (){
  if [ -z "${2}" ]; then
   echo "ERROR: missing argument for option '${1}'" 1>&2
   echo "see pantagruel --help for more details" 1>&2
   exit 1
  fi
}

ARGS=`getopt --options "h" --longoptions "ncpus:,mem:,wth:,hpctype:,chunksize:,parallelflags:,fwdenv:" --name "${hpcscript}" -- "$@"`

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
  export ncpus=${defncpus}
fi
if [ -z "${mem}" ] ; then
  export mem=${defmem}
fi
if [ -z "${wth}" ] ; then
  export wth=${defwth}
fi
if [ -z "${hpctype}" ] ; then
  export hpctype=${defhpctype}
fi
if [ -z "${chunksize}" ] ; then
  export chunksize=${defchunksize}
fi
if [ ! -z "${parallelflags}" ] ; then
  [ "${hpctype}" == 'PBS' ] && export parallelflags=":${parallelflags}"
  [ "${hpctype}" == 'LSF' ] && export parallelflags=" span[${parallelflags}]"
fi