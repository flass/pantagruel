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

checkptgversion
checkfoldersafe ${comparerecs}

###########################################
## 09. compare gene evolution scenarios
###########################################

mkdir -p ${compoutdir}/


## analyse of co-evolution!!
if [ -z $parsedreccolid ] ; then
  parsedreccolid=1
fi

## look for correlated gene lineage histories through identification of matching speciation and transfer events

### OPTION: exclude oldest species tree branches to avoid unspecific matches (and speed-up search):
if [ ! -z "${maxreftreeheight}" ] ; then
  step1="restricting events to those younger than age ${maxreftreeheight} on the species tree"
  echo ${step1}
  # e.g.: maxreftreeheight=0.5
  exclbrlist=${coretree}/branches_older_than_${maxreftreeheight}
  python2.7 ${ptgscripts}/list_branches.py --intree ${speciestreeBS}.lsd_internalPopulations.nwk --root_age 1.0 --older_than ${maxreftreeheight} --out ${exclbrlist}
  exclbr=$(cat ${exclbrlist} | tr '\n' ',' | sed -e "s/,$/\n/g" | sed -e "s/,/', '/g")
  
  if [ "${resumetask}" == 'true' ] ; then
    # varify if new tables were already built in which case revert to initial state
	if [ ! -z "$(sqlite3 ${sqldb} ".tables" | grep 'gene_lineage_events_full')" ] ; then
	  sqlite3 ${sqldb} """
	  DROP TABLE gene_lineage_events;
	  ALTER TABLE gene_lineage_events_full RENAME TO gene_lineage_events;
	  DROP INDEX IF EXISTS genelineev_rlocds;
      DROP INDEX IF EXISTS genelineev_rlocds_hash;
      DROP INDEX IF EXISTS genelineev_eventid;
      DROP INDEX IF EXISTS genelineev_eventid_hash;
      DROP INDEX IF EXISTS genelineev_freq;"""
    fi
  fi
  # create smaller table with only desired event
  sqlite3 ${sqldb} """
  ALTER TABLE gene_lineage_events RENAME TO gene_lineage_events_full;
  CREATE TABLE gene_lineage_events AS  
  (SELECT gene_lineage_events_full.* 
   FROM gene_lineage_events_full 
   INNER JOIN species_tree_events USING (event_id)
   INNER JOIN species_tree on rec_branch_id=branch_id
   WHERE branch_name NOT IN ('${exclbr}') )
   AND reconciliation_id=${parsedreccolid}
  ;
  CREATE INDEX genelineev_rlocds ON gene_lineage_events (replacement_label_or_cds_code);
  CREATE INDEX genelineev_rlocds_hash ON gene_lineage_events USING HASH (replacement_label_or_cds_code);
  CREATE INDEX genelineev_eventid ON gene_lineage_events (event_id);
  CREATE INDEX genelineev_eventid_hash ON gene_lineage_events USING HASH (event_id);
  CREATE INDEX genelineev_freq ON gene_lineage_events (freq);
  ALTER TABLE gene_lineage_events ADD PRIMARY KEY (replacement_label_or_cds_code, event_id);
  ANALYZE;
  .quit"""
  checkexec "failed ${step1}" "completed ${step1}\n"

fi
### end OPTION

# collect data
# BEWARE: GENERATES AN AWFUL LOT OF DATA, PREPARE DISK SPACE ACCORDINGLY
# indication: with defaults settings evtypematch='ST'; minevfreqmatch=0.5; minjoinevfreqmatch=1.0; maxreftreeheight=0.5; recsamplesize=1000
# on a 880 Enterobacteriaceae dataset, results in ~300 GB output (made to be split into ~1GB files)
python2.7 $ptgscripts/compare_collapsedALE_scenarios.py --events_from_postgresql_db ${sqldbname} --nrec_per_sample ${recsamplesize} \
 --event_type ${evtypematch} --min_freq ${minevfreqmatch} --min_joint_freq ${minjointevfreqmatch} --threads 8 \
 --dir_table_out ${compoutdir} &> ${ptglogs}/compare_collapsedALE_scenarios.${parsedreccol}.log &

if [ ! -z ${somehowusingpostgres} ] ; then
#### NOT IMPLEMENTED YET IN SQLite
# load data in database, adding mention of reconciliation_id to ensure events are not matched across collections
export assocoutdir=${compoutdir}/gene_lineage_assocations/between_fams_scores
${ptgscripts}/pantagruel_sqlitedb_load_coevolution_scores.py ${sqldb} ${assocoutdir} ${parsedreccolid}
####
fi