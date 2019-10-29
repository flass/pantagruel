#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019

if [ -z "$1" ] ; then
  echo "missing mandatory parameter: pantagruel config file"
  echo "Usage: $0 ptg_env_file [hpc_type:{PBS(default)|LSF}]"
  exit 1
fi
envsourcescript="$1"
if [ -z "$2" ] ; then 
  hpctype='PBS'
else
  hpctype="$2"
fi
envsourcescript="$1"
source ${envsourcescript}

############################
## 06.3 Bayesian gene trees
############################

## run mrbayes on collapsed alignments
export nexusaln4chains=${colalinexuscodedir}/${collapsecond}/collapsed_alns
export mboutputdir=${bayesgenetrees}/${collapsecond}
mkdir -p ${mboutputdir}
nchains=4
nruns=2
ncpus=$(( $nchains * $nruns ))
wth=24
tasklist=${nexusaln4chains}_ali_list
rm -f $tasklist
${ptgscripts}/lsfullpath.py "${nexusaln4chains}/*nex" > $tasklist
dtag=$(date +"%Y-%m-%d-%H-%M-%S")
append=''

if [ "${resumetask}" == 'true' ] ; then
  #~ # following lines are for resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
  for nfaln in $(cat $tasklist) ; do
   bnaln=$(basename $nfaln)
   bncontre=${bnaln/nex/mb.con.tre}
   if [ ! -e ${mboutputdir}/${bncontre} ] ; then
    echo $nfaln
   fi
  done > ${tasklist}_resumetask_${dtag}
  tasklist=${tasklist}_resumetask_${dtag}
  append='append=yes'
fi
qsubvar="tasklist=${tasklist}, outputdir=${mboutputdir}, mbmcmcpopt='Nruns=${nruns} Ngen=2000000 Nchains=${nchains}', mbmcmcopt='${append}'"

Njob=`wc -l ${tasklist} | cut -f1 -d' '`
chunksize=1000
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
for jobrange in ${jobranges[@]} ; do
 echo $jobrange $qsubvar
 dlogs=${ptglogs}/mrbayes/${chaintype}_mrbayes_trees_${collapsecond}_${dtag}_${jobrange}
 mkdir -p ${dlogs}/
 
 case "$hpctype" in
    'PBS') 
      subcmd="qsub -J ${jobrange} -N mb_panterodb -l select=1:ncpus=${ncpus}:mem=16gb -l walltime=${wth}:00:00 \
	  -o ${dlogs} -v \"$qsubvar\" ${ptgscripts}/mrbayes_array_PBS.qsub"
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
	  subcmd="bsub -J "raxml_gene_trees_$(basename $cdsfam2phylo)[$jobrange]" -q ${bqueue} \
	  -R \"select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]\" -n${ncpus} -M${memmb} \
	  -o ${dlogs}/mrbayes.%J.%I.o -e ${dlogs}/mrbayes.%J.%I.e -env \"$qsubvars\" ${ptgscripts}/mrbayes_array_PBS.bsub"
	  ;;
	*)
	  echo "Error: high-performance computer system '$hpctype' is not supported; exit now"
	  exit 1;;
    echo "${subcmd}"
    eval "${subcmd}"
  esac
done
