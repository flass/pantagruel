#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 15 Jan 2019

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file [ncpus=12] [mem=124gb] [walltimehours=24] [parallelflags='']" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

if [ ! -z "$2" ] ; then
  ncpus="$2"
else
  ncpus=12
fi
if [ ! -z "$3" ] ; then
  mem="$3"
else
  mem=124gb
fi
if [ ! -z "$4" ] ; then
  wth="$4"
else
  wth=24
fi
if [ ! -z "$5" ] ; then
  parallelflags=":$5"
  # for instance, parallelflags="mpiprocs=1:ompthreads=8"
  # will force the use of OpenMP multi-threading instead of default MPI parallelism
  # note this is not standard and might not be the default on your HPC system;
  # also the flags may be different depending on the HPC system config.
  # you should get in touch with your system admins to know the right flags
else
  parallelflags=""
fi



################################################################################
## 06.4 Convert format of Bayesian gene trees and replace species by populations
################################################################################

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A2: 
  echo "Error: not implemented yet:"
  echo "must generalize the script replace_species_by_pop_in_gene_trees.py so to only convert the format of gene tree chains, not replacing anything in them"
  # PBS-submitted parallel job
  #~ qsub -N replSpePopinGs -l select=1:ncpus=${ncpus}:mem=64gb,walltime=24:00:00 -o $repllogd -j oe -V << EOF
  #~ module load python
  #~ python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${tasklist} -o ${coltreechains} --threads=${ncpus} --reuse=0 --verbose=0 --logfile=${repllogs}_${replrun}.log &
  #~ EOF
  exit 1
  
  #### end OPTION A2: 
else
  #### OPTION B2: collapsed rake clades in gene trees need to be replaced by mock population leaves
  #### will edit collapsed gene trees to attribute an (ancestral) species identity to the leafs representing collapsed clades = pre-reconciliation of recent gene history
  if [ -z ${replacecolid} ] ; then
   export replacecolid=1
  fi
  mkdir -p ${coltreechains}/${collapsecond}

  ## edit the gene trees, producing the corresponding (potentially collapsed) species tree based on the 'time'-tree backbone
  mkdir -p ${ptgdb}/logs/replspebypop
  export tasklist=${bayesgenetrees}_${collapsecond}_mbrun1t_list
  ls ${bayesgenetrees}/${collapsecond}/*run1.t > $tasklist
  repllogd=${ptgdb}/logs/replspebypop
  repllogs=${repllogd}/replace_species_by_pop_in_gene_trees
  replacecoldate=$(date +%Y-%m-%d)
  echo -e "${replacecolid}\t${replacecoldate}" > ${genetrees}/replacecol
  
  export dtag="$(date +'%Y%m%d-%H%M%S')"
  if [ "${resumetask}" == 'true' ] ; then
    #~ # following lines are for resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
    for nfrun1t in $(cat $tasklist) ; do
     bnrun1t=$(basename $nfrun1t)
     bnGtre=${bnrun1t/.mb.run1.t/-Gtrees.nwk}
     bnStre=${bnrun1t/.mb.run1.t/-Stree.nwk}
     if [[ ! -e ${coltreechains}/${collapsecond}/${colmethod}/${bnGtre} || ! -e ${coltreechains}/${collapsecond}/${colmethod}/${bnStre} ]] ; then
       echo ${nfrun1t}
     fi
    done > ${tasklist}_resumetasklist
    export tasklist=${tasklist}_resumetasklist
  fi
  # PBS-submitted parallel job
  # divide run in small chunks o be run in different jobs
  chunksize=100
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
  qsub${arrayspec} -N replSpePopinGs -l select=1:ncpus=${ncpus}:mem=${mem}${parallelflags},walltime=${wth}:00:00 -o ${repllogd} -j oe -V ${ptgscripts}/replSpePopinGs_array_PBS.qsub
  
  ## load these information into the database
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh ${database} ${sqldb} ${colalinexuscodedir} ${coltreechains} ${collapsecond} ${colmethod} "${collapsecriteriondef}" ${collapsecolid} ${replacecolid} ${collapsecoldate} ${replacecoldate}

fi
#### end OPTION B2
