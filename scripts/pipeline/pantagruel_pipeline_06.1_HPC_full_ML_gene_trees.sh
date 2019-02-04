#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019


if [ -z "$2" ] ; then echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgroot="$1"    # source folder where to create the database
export ptgdbname="$2"  # database anme (will notably be the name of the top folder)
export ptgdb=${ptgroot}/${ptgdbname}

envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
source ${envsourcescript}


#############################################################
## 06.1 Full ML gene trees on HPC
#############################################################

## compute first pass of gene trees with RAxML, using rapid bootstrap to estimate branch supports
allcdsfam2phylo=($(ls ${protali}/cdsfams_*))
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
qsub -J $jobrange -l walltime=${wt}:00:00 -l select=1:ncpus=${ncpus}:mem=${mem}gb -N raxml_gene_trees_$(basename $cdsfam2phylo) -o $ptglogs/raxml/gene_trees -j oe -v "$qsubvars" ${ptgscripts}/raxml_array_PBS.qsub
done
done

