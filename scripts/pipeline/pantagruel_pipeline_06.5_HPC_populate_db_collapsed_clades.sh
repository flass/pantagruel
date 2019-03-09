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
## 06.5 Populate database withh all gene tree results
################################################################################

if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  #### OPTION A2: pretty much nothing to do
  echo "Error: not implemented yet exit now"
  exit 1
  
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
