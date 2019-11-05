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

ARGS=`getopt --options "h" --longoptions "ncpus:,mem:,wth:,hpctype:,chunksize:,parallelflags:,fwdenv:" --name "pantagruel_pipeline_06.2_HPC_collapse_gene_trees.sh" -- "$@"`

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
  mem=32
fi
if [ -z "${wth}" ] ; then
  wth=4
fi
if [ -z "${hpctype}" ] ; then
  hpctype='PBS'
fi
if [ -z "${chunksize}" ] ; then
  chunksize=3000
fi
if [ ! -z "${parallelflags}" ] ; then
  [ "${hpctype}" == 'PBS' ] && parallelflags=":${parallelflags}"
  [ "${hpctype}" == 'LSF' ] && parallelflags=" span[${parallelflags}]"
fi
#if [ ! -z "${fwdenv}" ] ; then
#  [ "${hpctype}" == 'PBS' ] && fwdenv="-v '${fwdenv}'"
#  [ "${hpctype}" == 'LSF' ] && fwdenv="-env '${fwdenv}'"
#fi

#############################################################
## 06.2 Gene tree collapsing on HPC
#############################################################

# first collate the list of families for which no gene tree was computed 
# due to containing too few non-identical sequences 
# and will be discarded for the remamining analyses
smallfams=${genetrees}/reduced_alignment_is_too_small_fams
touch ${smallfams}
cat ${genetrees}/bulk/*.smallreducedali > ${smallfams}

echo "$(wc -l ${smallfams}) gene families were found to contain to few non-identical CDSs for a gene tree to be computed"

#### OPTION: edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  # no collapsing, just convert the alignments from fasta to nexus
  for aln in `ls ${alifastacodedir}` ; do
    convalign -i fasta -e nex -t dna nexus ${alifastacodedir}/$alnfa 
  done
  mv ${alifastacodedir}/*nex ${colalinexuscodedir}/
  
else
 if [[ "${chaintype}" == 'collapsed' ]] ; then
  if [ -z ${collapsecolid} ] ; then
    collapsecolid=1
  fi
  if [[ "$collapseCladeParams" != 'default' ]] ; then
    eval "$collapseCladeParams"
    # e.g.:  eval 'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'
  fi
  
  ## detect clades to be collapsed in gene trees
  collapsecriteriondef="--clade_stem_conds=[('$criterion','>=',$cladesupp)] --within_clade_conds=[('$withinfun','$criterion','<=',$subcladesupp,-1),('max','$criterion','<',$cladesupp,-1)]"
  mkdir -p ${colalinexuscodedir}/${collapsecond}
  echo "${collapsecriteriondef}" > ${colalinexuscodedir}/${collapsecond}.collapse_criterion_def
  mlgenetreelist=${mlgenetrees%*s}_list
  ${ptgscripts}/lsfullpath.py "${mlgenetrees}/${mainresulttag}/*" > ${mlgenetreelist}
  
  # accomodate with possible upper limit on number of tasks in an array job; assume chunks of 3000 tasks are fine
  Njob=`wc -l ${mlgenetreelist} | cut -f1 -d' '`
  jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))

  # In the job submission commands below, some lines are specific to the HPC system 
  # and environmnt on which the script was developped:
  #   module load anaconda3/personal
  #   source activate env_python2
  # In other environments, other methods may be used to access the required Python packages.
  # To emulate this on other systems, it is best to use anaconda to build your own environment
  # where you make sure you are using python 2.7 and have all required packages installed with it.
  # You can do so using the following command:
  # conda create -n env_python2 python=2.7 python-igraph biopython bcbio-gff scipy

  for jobrange in ${jobranges[@]} ; do
    beg=`echo $jobrange | cut -d'-' -f1`
    tail -n +${beg} ${mlgenetreelist} | head -n ${chunksize} > ${mlgenetreelist}_${jobrange}
	
    qsubvars="mlgenetreelist=${mlgenetreelist}_${jobrange},ptgscripts,cdsalifastacodedir,ncpus,colalinexuscodedir,collapsecond,mlgenetrees"
    if [ ! -z "${fwdenv}" ] ; then
	  qsubvars="${qsubvars},${fwdenv}"
	fi
	
	case "${hpctype}" in
      'PBS')
         subcmd="qsub -N 'mark_unresolved_clades' -l select=1:ncpus=${ncpus}:mem=${mem}gb,walltime=${wth}:00:00 \
          -o ${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.log -j oe \
		  -v \"${qsubvars}\" ${ptgscripts}/mark_unresolved_clades.qsub"
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
	    nflog="${ptglogs}/mark_unresolved_clades.${collapsecond}_${jobrange}.%J.o"
	    subcmd="bsub -J 'mark_unresolved_clades' -q ${bqueue} \
	     -R \"select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]\" \
         -n${ncpus} -M${memmb} -env \"$qsubvars\" \
         -o ${nflog} -e ${nflog} ${ptgscripts}/mark_unresolved_clades.bsub"
	    ;;
	  *)
	    echo "Error: high-performance computer system '$hpctype' is not supported; exit now"
	    exit 1;;
    esac
	echo "${subcmd}"
	eval "${subcmd}"
  done
  export collapsecoldate=$(date +%Y-%m-%d)
  echo -e "${collapsecolid}\t${collapsecoldate}" > ${genetrees}/collapsecol
  
 else
  echo "Error: incorrect value for variable chaintype: '${chaintype}'"
  exit 1
 fi
fi
#### end OPTION
