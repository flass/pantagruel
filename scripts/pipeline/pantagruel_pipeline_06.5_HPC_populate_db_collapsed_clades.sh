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
## 06.5 Populate database withh all gene tree results
################################################################################

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A2: 
  echo "Error: not implemented yet:"
  echo "nothing to do; exit now"
  exit 0
  
  #### end OPTION A2: 
else
  #### OPTION B2: rake clades in gene trees were collapsed and later replaced by mock population leaves
  #### will feed data relative to these operation to the SQL databse

  ## edit the gene trees, producing the corresponding (potentially collapsed) species tree based on the 'time'-tree backbone
  collapsecolid=$(cut -f1 ${genetrees}/collapsecol)
  collapsecoldate=$(cut -f2 ${genetrees}/collapsecol)
  collapsecriteriondef="$(cut -f3 ${genetrees}/collapsecol)"
  replacecolid=$(cut -f1 ${genetrees}/replacecol)
  replacecoldate=$(cut -f2 ${genetrees}/replacecol)
    
  ## load these information into the database
  ${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_collapsed_clades.sh ${database} ${sqldb} ${colalinexuscodedir} ${coltreechains} ${collapsecond} ${colmethod} "${collapsecriteriondef}" ${collapsecolid} ${replacecolid} ${collapsecoldate} ${replacecoldate}

fi
#### end OPTION B2
