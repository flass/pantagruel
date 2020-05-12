#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file [gene_fam_list]" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

if [ ! -z "$2" ] ; then
  export genefamlist="$2"
fi

###############################################
## 07. Gene tree / Species tree reconciliations
###############################################

echo "Will use the reconciliation method: ${recmethod}"
if [ "${recmethod}" != 'GeneRax' ] ; then
  echo "This is the 'usingGeneRax' branch of Pantagruel pipeline; only 'GeneRax' reconciliation method is supported."
  echo "You can edit the environment file '${envsourcescript}' to set the value of the variable with the command \`export recmethod='GeneRax'\`;"
  echo "or you can use the 'master' branch of the pipeline for other reconciliation methods."
  echo "Pantagruel will exit now."
  exit 1
fi
if [ "${recmethod}" == 'ALE' ] ; then
  ${ptgscripts}/pipeline/pantagruel_pipeline_07_ALE_reconciliations.sh ${envsourcescript}
elif [ "${recmethod}" == 'ecceTERA' ] ; then
  ${ptgscripts}/pipeline/pantagruel_pipeline_07_ecceTERA_reconciliations.sh ${envsourcescript}
  elif [ "${recmethod}" == 'GeneRax' ] ; then
  ${ptgscripts}/pipeline/pantagruel_pipeline_07_GeneRax_reconciliations.sh ${envsourcescript}
fi
checkexec "Could not complete reconciliations with ${recmethod}" "Successfully reconciled gene trees with ${recmethod}"
