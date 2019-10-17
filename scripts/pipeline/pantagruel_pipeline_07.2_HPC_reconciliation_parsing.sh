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

# safer specifying the reconciliation collaction to parse; e.g.: 'ale_collapsed_dated_1'
if [ ! -z "$2" ] ; then
  reccol="$2"
else
  # if not inferred from the record of the last reconciliation computation
  reccol=$(cut -f4 ${alerec}/reccol)
fi

if [ ! -z "$3" ] ; then
  ncpus="$3"
else
  ncpus=8
fi
if [ ! -z "$4" ] ; then
  mem="$"
else
  mem=124gb
fi
if [ ! -z "$5" ] ; then
  wth="$5"
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

outrecdir=${recs}/${collapsecond}/${replmethod}/${reccol}

mkdir -p ${parsedrecs}
reclist=$outrecdir/${reccol}_${rectype}_uml_rec_list
${ptgscripts}/lsfullpath.py "${outrecdir}/${reccol}_${rectype}/*ml_rec" > $reclist
 
## normalise the species tree branch labels across gene families
## and look for correlated transfer events across gene families

# PBS-submitted parallel job
parsecollogd=${ptgdb}/logs/parsecol
parsecollogs=${parsecollogd}/parse_collapsedALE_scenarios.og

# In the job submission commands below, some lines are specific to the HPC system 
# and environmnt on which the script was developped:
#   module load anaconda3/personal
#   source activate env_python2
# In other environments, other methods may be used to access the required Python packages.
# To emulate this on other systems, it is best to use anaconda to build your own environment
# where you make sure you are using python 2.7 and have all required packages installed with it.
# You can do so using the following command:
# conda create -n env_python2 python=2.7 python-igraph biopython bcbio-gff scipy
ALEsourcenote=$(cut -f3 ${alerec}/reccol)

#### NOTE
## here for simplicity only the log variables refering to parsed reconciliations (parsedreccol, parsedreccolid, parsedreccoldate) are recorded in the database
## but not the log variables refering to the actual reconciliations (reccol, reccolid, reccoldate)
####

qsub -N parseColALE -l select=1:ncpus=${ncpus}:mem=${mem},walltime=${wth}:00:00 -o ${parsecollogs} -j oe -V << EOF
module load anaconda3/personal
source activate env_python2
python2.7 ${ptgscripts}/parse_collapsedALE_scenarios.py --rec_sample_list ${reclist} \
 --populations ${speciestree/.full/}_populations --reftree ${speciestree}.lsd.nwk \
 --dir_table_out ${parsedrecs} --evtype ${evtypeparse} --minfreq ${minevfreqparse} \
 --threads 8

export parsedreccoldate=$(date +%Y-%m-%d)
echo -e "${parsedreccolid}\t${parsedreccoldate}" > ${alerec}/parsedreccol

## store reconciliation parameters and load parsed reconciliation data into database
${ptgscripts}/pantagruel_sqlitedb_phylogeny_populate_reconciliations.sh "${database}" "${sqldb}" "${parsedrecs}" "${ALEversion}" "${ALEalgo}" "${ALEsourcenote}" "${parsedreccol}" "${parsedreccolid}" "${parsedreccoldate}"
EOF
