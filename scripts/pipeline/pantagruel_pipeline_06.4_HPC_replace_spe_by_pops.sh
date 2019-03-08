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
   replacecolid=1
  fi
  mkdir -p ${coltreechains}/${collapsecond}

  ## edit the gene trees, producing the corresponding (potentially collapsed) species tree based on the 'time'-tree backbone
  mkdir -p ${ptgdb}/logs/replspebypop
  tasklist=${bayesgenetrees}_${collapsecond}_mbrun1t_list
  ls ${bayesgenetrees}/${collapsecond}/*run1.t > $tasklist
  repllogd=${ptgdb}/logs/replspebypop
  repllogs=${repllogd}/replace_species_by_pop_in_gene_trees
  replrun=$(date +'%d%m%Y')

  if [ "${resumetask}" == 'true' ] ; then
    #~ # following lines are for resuming after a stop in batch computing, or to collect those jobs that crashed (and may need to be re-ran with more mem/time allowance)
    for nfrun1t in $(cat $tasklist) ; do
     bnrun1t=$(basename $nfrun1t)
     bnGtre=${bnrun1t/mb.run1.t/Gtrees.nwk}
     bnStre=${bnrun1t/mb.run1.t/Stree.nwk}
     if [[ ! -e ${coltreechains}/${collapsecond}/${colmethod}/${bnGtre} || ! -e ${coltreechains}/${collapsecond}/${colmethod}/${bnStre} ]] ; then
       echo ${nfrun1t}
     fi
    done > ${tasklist}_resumetask_${dtag}
    tasklist=${tasklist}_resumetask_${dtag}
  fi
  # PBS-submitted parallel job
  ncpus=8
  qsub -N replSpePopinGs -l select=1:ncpus=${ncpus}:mem=96gb,walltime=72:00:00 -o ${repllogd} -j oe -V << EOF
  module load python
  python ${ptgscripts}/replace_species_by_pop_in_gene_trees.py -G ${tasklist} -c ${colalinexuscodedir}/${collapsecond} -S ${speciestree}.lsd.newick -o ${coltreechains}/${collapsecond} \
   --populations=${speciestree%.*}_populations --population_tree=${speciestree%.*}_collapsedPopulations.nwk --population_node_distance=${speciestree%.*}_interNodeDistPopulations \
   --dir_full_gene_trees=${mlgenetrees}/rootedTree --method=${colmethod} --threads=${ncpus} --reuse=0 --max.recursion.limit=12000 --logfile=${repllogs}_${replrun}.log
  if [ $? != 0 ] ; then
    echo "ERROR: ${ptgscripts}/replace_species_by_pop_in_gene_trees.py call returned wth an error status ; quit now" 1>&2
    exit 1
  fi
  export replacecoldate=$(date +%Y-%m-%d)
  echo -e "${replacecolid}\t${replacecoldate}" > ${genetrees}/replacecol
  ## load these information into the database
  collapsecolid=$(cut -f1 ${genetrees}/collapsecol)
  collapsecoldate=$(cut -f2 ${genetrees}/collapsecol)
  collapsecriteriondef="$(cut -f3 ${genetrees}/collapsecol)"
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh ${database} ${sqldb} ${colalinexuscodedir} ${coltreechains} ${collapsecond} ${colmethod} "${collapsecriteriondef}" ${collapsecolid} ${replacecolid} ${collapsecoldate} ${replacecoldate}
EOF


fi
#### end OPTION B2
