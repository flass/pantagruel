#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019


if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}


#############################################################
## 06.1 Full ML gene trees on HPC
#############################################################

mkdir -p ${ptglogs}/raxml/gene_trees/ ${$mlgenetrees}/

## compute first pass of gene trees with RAxML, using rapid bootstrap to estimate branch supports
rm -f ${protali}/cdsfams_*_aln_list
sortminsizes=($(for famset in `ls ${protali}/*cdsfams_*` ; do echo $famset | sed -e 's/.*cdsfams_minsize\([0-9]\+.*\)/\1/' ; done | sort -n))
allcdsfam2phylo=($(for famminsize in ${sortminsizes[@]} ; do ls ${protali}/*${famminsize} ; done))
allmems=(4 8 32 64)
allwalltimes=(24 24 72 72)
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
qsubvars="tasklist=$tasklist,outputdir=$mlgenetrees,reducedaln=true,nbthreads=${ncpus}"
Njob=`wc -l ${tasklist} | cut -f1 -d' '`
# accomodate with possible upper limit on number of tasks in an array job; assume chunks of 3000 tasks are fine
chunksize=3000
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
for jobrange in ${jobranges[@]} ; do
echo "jobrange=$jobrange"
qsub -J $jobrange -l walltime=${wt}:00:00 -l select=1:ncpus=${ncpus}:mem=${mem}gb -N raxml_gene_trees_$(basename $cdsfam2phylo) -o ${ptglogs}/raxml/gene_trees -j oe -v "$qsubvars" ${ptgscripts}/raxml_array_PBS.qsub
done
done

