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
source ${envsourcescript}


#############################################################
## 06.1 Full ML gene trees on HPC
#############################################################

mkdir -p ${ptglogs}/raxml/gene_trees/ ${$mlgenetrees}/

## compute first pass of gene trees with RAxML, using rapid bootstrap to estimate branch supports
if [[ "${resumetask}" == "true" && ! -z "$(ls ${genetrees}/cdsfams_*_aln_list 2> /dev/null)" ]] ; then
  echo "will re-use previously established lists of gene families for which to compute trees" 
  allcdsfam2phylo=($(ls ${genetrees}/cdsfams_*_aln_list))
else
  rm -f ${genetrees}/cdsfams_*_aln_list
  sortminsizes=($(for famset in `ls ${genetrees}/*cdsfams_*` ; do echo $famset | sed -e 's/.*cdsfams_minsize\([0-9]\+.*\)/\1/' ; done | sort -n))
  allcdsfam2phylo=($(for famminsize in ${sortminsizes[@]} ; do ls ${genetrees}/*${famminsize} ; done))
fi
  
allmems=(4 8 32 64)
allwalltimes=(12 24 48 72)
allncpus=(4 4 4 16)
for i in ${!allcdsfam2phylo[@]} ; do
cdsfam2phylo=${allcdsfam2phylo[$i]} ; mem=${allmems[$i]} ; wt=${allwalltimes[$i]} ; ncpus=${allncpus[$i]}
echo "cdsfam2phylo=${cdsfam2phylo} ; mem_per_core=${mem}gb ; walltime=${wt}:00:00 ; ncpus=${ncpus}"
tasklist=${cdsfam2phylo}_aln_list
rm -f $tasklist ; for fam in `cut -f1 $cdsfam2phylo` ; do ls ${cdsalifastacodedir}/${fam}.codes.aln >> $tasklist ; done
if [ "$(wc -l $cdsfam2phylo | cut -f1 -d' ')" -lt "$(wc -l $tasklist | cut -f1 -d' ')" ] ; then 
  >&2 echo "ERROR $(dateprompt): missing gene family alignments; please fix the list '$tasklist' or the content of folder '$alifastacodedir/' before continuing."
  exit 1
fi
rm -f ${tasklist}* ; for fam in $(cut -f1 ${cdsfam2phylo}) ; do
  aln=${cdsalifastacodedir}/${fam}.codes.aln
  if [[ "${resumetask}" == "true" ]] ; then
    if [[ ! -s ${mlgenetrees}/${mainresulttag}/RAxML_${mainresulttag}.${fam}.codes ]] ; then
      ls ${aln} >> ${tasklist}_resume
    fi
  fi
  ls ${aln} >> ${tasklist}
done
if [[ "${resumetask}" == "true" ]] ; then
  if [[ -s ${tasklist}_resume ]] ; then
    echo "Resume task 6: $(wc -l ${tasklist}_resume | cut -d' ' -f1) ML trees left to infer"
  else
    echo "Resume task 6: all ML tree built; skip ML tree building"
  fi
  tasklist=${tasklist}_resume
fi



qsubvars="tasklist=$tasklist,outputdir=$mlgenetrees,reducedaln=true,nbthreads=${ncpus}"
Njob=`wc -l ${tasklist} | cut -f1 -d' '`
# accomodate with possible upper limit on number of tasks in an array job; assume chunks of 3000 tasks are fine
chunksize=1000
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
for jobrange in ${jobranges[@]} ; do
echo "jobrange=$jobrange"
  case "$hpctype" in
    'PBS') 
      qsub -J $jobrange -l walltime=${wt}:00:00 -l select=1:ncpus=${ncpus}:mem=${mem}gb -N raxml_gene_trees_$(basename $cdsfam2phylo) \
	  -o ${ptglogs}/raxml/gene_trees -j oe -v "$qsubvars" ${ptgscripts}/raxml_array_PBS.qsub
	  ;;
	'LSF')
	  if [ ${wt} -le 12 ] ; then
	    bqueue='normal'
	  elif [ ${wt} -le 48 ] ; then
	    bqueue='long'
	  else
	    bqueue='basement'
	  fi
	  memmb=$((${mem} * 1000)) 
	  bsub -J "raxml_gene_trees_$(basename $cdsfam2phylo)[$jobrange]" -q ${bqueue} -R "select[mem>${memmb}] rusage[mem=${memmb}] span[hosts=1]" \
	  -n${ncpus} -M${memmb} -o ${ptglogs}/raxml/gene_trees/raxml_gene_trees_$(basename $cdsfam2phylo).%J.%I.o \
	  -e ${ptglogs}/raxml/gene_trees/raxml_gene_trees_$(basename $cdsfam2phylo).%J.%I.e -env "$qsubvars" ${ptgscripts}/raxml_array_LSF.bsub
	  ;;
	*)
	  echo "Error: high-performance computer system '$hpctype' is not supported; exit now"
	  exit 1;;
  esac
done
done

