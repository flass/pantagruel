#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file [ncpus=8] [mem=124gb] [walltimehours=24]" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

if [ ! -z "$2" ] ; then
  ncpus="$2"
else
  ncpus=8
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

######################################################
## 07.2 Parse gene tree / Species tree reconciliations
######################################################

### parse the inferred scenarios
# parameters to be set
if [ -z $parsedreccolid ] ; then
  parsedreccolid=1
fi
# derived parameters
export parsedreccol=${reccol}_parsed_${parsedreccolid}
export parsedrecs=${alerec}/parsed_recs/${parsedreccol}

mkdir -p ${parsedrecs}
reclist=$outrecdir/ale_collapsed_${rectype}_uml_rec_list
${ptgscripts}/lsfullpath.py "${outrecdir}/ale_collapsed_${rectype}/*ml_rec" > $reclist
 
## normalise the species tree branch labels across gene families
## and look for correlated transfer events across gene families

# PBS-submitted parallel job
parsecollogd=${ptgdb}/logs/parsecol
parsecollogs=${parsecollogd}/parse_collapsedALE_scenarios.og

qsub -N parseColALE -l select=1:ncpus=${ncpus}:mem=${mem},walltime=${wth}:00:00 -o ${parsecollogs} -j oe -V << EOF
module load python
python ${ptgscripts}/parse_collapsedALE_scenarios.py --rec_sample_list ${reclist} \
 --populations ${speciestree/.full/}_populations --reftree ${speciestree}.lsd.newick \
 --dir_table_out ${parsedrecs} --evtype ${evtypeparse} --minfreq ${minevfreqparse} \
 --threads 8

export parsedreccoldate=$(date +%Y-%m-%d)
echo -e "${parsedreccolid}\t${parsedreccoldate}" > ${alerec}/parsedreccol

## store reconciliation parameters and load parsed reconciliation data into database
${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_reconciliations.sh ${database} ${sqldb} ${parsedrecs} ${ALEversion} ${ALEalgo} ${ALEsourcenote} ${parsedreccol} ${parsedreccolid} ${parsedreccoldate}
EOF
