#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019

if [ -z "$2" ] ; then echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
export ptgdb=${ptgroot}/${ptgdbname}

envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
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
tasklist=${nexusaln4chains}_ali_list
rm -f $tasklist
${ptgscripts}/lsfullpath.py ${nexusaln4chains} > $tasklist

#~ # following lines are for resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
#~ alreadytrees=${mboutputdir}_list
#~ ${ptgscripts}/lsfullpath.py ${mboutputdir}/*con.tre > $alreadytrees
#~ alreadytasklist=${nexusaln4chains}_ali_list_done
#~ sed -e "s#${mboutputdir}/\(.\+\)\.mb\.con\.tre#${nexusaln4chains}/\1\.nex#g" $alreadytrees > $alreadytasklist
#~ sort $tasklist > $tasklist.sort
#~ sort $alreadytasklist > $alreadytasklist.sort
#~ dtag=$(date +"%Y-%m-%d-%H-%M-%S")
#~ comm -2 -3  $tasklist.sort $alreadytasklist.sort > ${tasklist}_todo_${dtag}
#~ Njob=`wc -l ${tasklist}_todo_${dtag} | cut -f1 -d' '`
#~ qsubvar="mbversion=3.2.6, tasklist=${tasklist}_todo_${dtag}, outputdir=${mboutputdir}, mbmcmcpopt='Nruns=${nruns} Ngen=2000000 Nchains=${nchains}'"
# otherwise could just use $tasklist
dtag=$(date +"%Y-%m-%d-%H-%M-%S")
Njob=`wc -l ${tasklist} | cut -f1 -d' '`
chunksize=1000
jobranges=($(${ptgscripts}/get_jobranges.py $chunksize $Njob))
qsubvar="mbversion=3.2.6, tasklist=${tasklist}, outputdir=${mboutputdir}, mbmcmcpopt='Nruns=${nruns} Ngen=2000000 Nchains=${nchains}'"
for jobrange in ${jobranges[@]} ; do
 echo $jobrange $qsubvar
 qsub -J $jobrange -N mb_panterodb -l select=1:ncpus=${ncpus}:mem=16gb -o ${ptglogs}/mrbayes/collapsed_mrbayes_trees_${dtag}_${jobrange} -v "$qsubvar" ${ptgscripts}/mrbayes_array_PBS.qsub
done
