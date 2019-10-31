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
 For running the script:\n
\t$0 ptg_config_file [ncpus=12] [mem=124] [wth=24] [hpctype='PBS'] [chunksize=100] [parallelflags='']\n
or for this help mesage:\n
\t$0 -h\n
\n
\tOption details:\n
\t ncpus:\tnumber of parallel threads.\n
\t mem:\tmax memory allowance, in gigabytes (GB) (default: 124).\n
\t wth:\tmax walltime, in hours (default: 24).\n
\t hpctype:\t'PBS' for Torque scheduling system (default), or 'LSF' for IBM schedulling system.\n
\t chunksize:\tnumber of families to be processed in one submitted job (default: 100)\n
\t parallelflags:\tadditional flags to be specified for worker nodes' resource requests\n
\t\t\t for instance, specifying parallelflags=\"mpiprocs=1:ompthreads=12\" (with hpctype='PBS')\n
\t\t\t will force the use of OpenMP multi-threading instead of default MPI parallelism.\n
\t\t\t Note that parameter declaration syntax is not standard and differ depending on your HPC system;\n
\t\t\t also the relevant parameter flags may be different depending on the HPC system config.\n
\t\t\t Get in touch with your system admins to know the relevant flags and syntax.\n
"

if [ "$1" == '-h' ] ; then
 echo -e ${usagemsg}
 exit 0
elif [ -z "$1" ] ; then
 echo -e "Error: missing mandatory parameter: pantagruel config file\n"
 echo -e ${usagemsg}
 exit 1
fi
envsourcescript="$1"
source ${envsourcescript}

if [ ! -z "$2" ] ; then
  export ncpus="$2"
else
  export ncpus=12
fi
if [ ! -z "$3" ] ; then
  mem="$3"
else
  mem=124
fi
if [ ! -z "$4" ] ; then
  wth="$4"
else
  wth=24
fi
if [ ! -z "$5" ] ; then
  hpctype="$5"
else
  hpctype='PBS'
fi
if [ ! -z "$6" ] ; then
  chunksize="$6"
else
  chunksize=100
fi
if [ ! -z "$7" ] ; then
  [ "$hpctype" == 'PBS' ] && parallelflags=":$7"
  [ "$hpctype" == 'LSF' ] && parallelflags="span[$7]"
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
fi
