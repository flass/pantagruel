#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$1" ] ; then echo "missing mandatory parameter: pantagruel config file" ; echo "Usage: $0 ptg_env_file" ; exit 1 ; fi
envsourcescript="$1"
source ${envsourcescript}

checkfoldersafe (){
  if [ -d ${1} ] ; then
    if [ ${runmode} != 'force' ] ; then
      echo "Task folder ${1} already exists; will stop here rther then overwritting data."
      echo "If you want previous data to be reased, use FORCE mode with -F|--FORCE option."
      exit 2
    else
      echo "Task folder ${1} already exists; FORCE mode is on: WILL ERASE the folder and write new result in its place"
      rm -r ${1}
    fi
  else
    mkdir ${1}
  fi
}

checkfoldersafe ${comparerecs}

###########################################
## 09. compare gene evolution scenarios
###########################################

mkdir -p ${compoutdir}/


## analyse of co-evolution!!
export minevfreqmatch=0.5
export minjoinevfreqmatch=1.0
if [ -z $parsedreccolid ] ; then
  parsedreccolid=1
fi

## look for correlated gene lineage histories through identification of matching speciation and transfer events

### OPTION: exclude oldest species tree branches to avoid unspecific matches (and speed-up search):
if [ ! -z maxreftreeheight ]
  # e.g.: maxreftreeheight=0.25
  exclbrlist=${coretree}/branches_older_than_${maxreftreeheight}
  python ${ptgscripts}/list_branches.py --intree ${speciestreeBS}.lsd_internalPopulations.nwk --root_age 1.0 --older_than ${maxreftreeheight} --out ${exclbrlist}
  exclbr=$(cat ${exclbrlist} | tr '\n' ',' | sed -e "s/,$/\n/g" | sed -e "s/,/', '/g")

  # create smaller table with only desired event
  sqlite3 ${sqldb} << EOF
  ALTER TABLE gene_lineage_events RENAME TO gene_lineage_events_full;
  CREATE TABLE gene_lineage_events AS  
  (SELECT gene_lineage_events_full.* 
   FROM gene_lineage_events_full 
   INNER JOIN species_tree_events USING (event_id)
   INNER JOIN species_tree on rec_branch_id=branch_id
   WHERE branch_name NOT IN ('${exclbr}') )
   AND reconciliation_id=${parsedreccolid}
  ;
  CREATE INDEX ON gene_lineage_events (replacement_label_or_cds_code);
  CREATE INDEX ON gene_lineage_events USING HASH (replacement_label_or_cds_code);
  CREATE INDEX ON gene_lineage_events (event_id);
  CREATE INDEX ON gene_lineage_events USING HASH (event_id);
  CREATE INDEX ON gene_lineage_events (freq);
  ALTER TABLE gene_lineage_events ADD PRIMARY KEY (replacement_label_or_cds_code, event_id);
  ANALYZE;
  .quit
EOF

fi
### end OPTION

# collect data
# BEWARE: GENERATES AN AWFUL LOT OF DATA, PREPARE DISK SPACE ACCORDINGLY
# indication: with defaults settings evtypeparse='ST'; minevfreqmatch=0.5; minjoinevfreqmatch=1.0; maxreftreeheight=0.25
# on a 880 Enterobacteriaceae dataset, results in ~300 GB output (made to be split into ~1GB files)
python $ptgscripts/compare_collapsedALE_scenarios.py --events_from_postgresql_db ${sqldbname} \
 --event_type ${evtypeparse} --min_freq ${minevfreqmatch} --min_joint_freq ${minjointevfreqmatch} --threads 8 \
 --dir_table_out ${compoutdir} &> $entlogs/compare_collapsedALE_scenarios.${parsedreccol}.log &

#### NOT IMPLEMENTED YET IN SQLite
# load data in database, adding mention of reconciliation_id to ensure events are not matched across collections
export assocoutdir=${compoutdir}/gene_lineage_assocations/between_fams_scores
$ptgscripts/pantagruel_sqlitedb_load_coevolution_scores.py ${sqldb} ${assocoutdir} ${parsedreccolid}
####
