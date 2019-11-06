#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019

# environment variables to be passed on to interface script:
export hpcscript=$(basename ${0})
hpcscriptdir=$(dirname ${0})
export defncpus=8
export defmem=96
export defwth=24
export defhpctype='PBS'
export defchunksize=100
export withpython='true'

source ${hpcscriptdir}/pantagruel_pipeline_HPC_common_interface.sh "${@}"

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

