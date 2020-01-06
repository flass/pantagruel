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

# safer specifying the parsed reconciliation collection to consider; e.g.: 'ale_collapsed_dated_1_parsed_1'
if [ -z "${reccol}" ] ; then
  # if not inferred from the record of the last reconciliation computation
  parsedreccol=$(cut -f3 ${alerec}/parsedreccol)
  [ -z "${reccol}" ] && echo "Error: cannot find parsed reconciliation collection as env variable \$parsedreccol is empty; exit now" && exit 1
fi

## look for correlated gene lineage histories through identification of matching speciation and transfer events

### OPTION: exclude oldest species tree branches to avoid unspecific matches (and speed-up search):
if [ ! -z "${maxreftreeheight}" ] ; then
  step1="restricting events to those younger than age ${maxreftreeheight} on the species tree"
  echo ${step1}
  # e.g.: maxreftreeheight=0.5
  exclbrlist=${coretree}/branches_older_than_${maxreftreeheight}
  python2.7 ${ptgscripts}/list_branches.py --intree ${speciestree}.lsd_internalPopulations.nwk --root_age 1.0 --older_than ${maxreftreeheight} --out ${exclbrlist}
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

step2="computing co-evolution scores"
# collect data
spacewarning="""BEWARE: GENERATES AN AWFUL LOT OF DATA, PREPARE DISK SPACE ACCORDINGLY
# indication: with defaults settings evtypematch='ST'; minevfreqmatch=0.5; minjoinevfreqmatch=1.0; maxreftreeheight=0.25; recsamplesize=1000
# on a 880 Enterobacteriaceae dataset, results in ~300 GB output (made to be split into ~1GB files)"""
echo ${step2}
echo ${spacewarning}
python2.7 $ptgscripts/compare_collapsedALE_scenarios.py --events_from_sqlite_db ${sqldb} --nrec_per_sample ${recsamplesize} \
 --event_type ${evtypematch} --min_freq ${minevfreqmatch} --min_joint_freq ${minjointevfreqmatch} --threads 8 \
 --dir_table_out ${compoutdir} &> ${ptglogs}/compare_collapsedALE_scenarios.${parsedreccol}.log
checkexec "failed ${step2}" "completed ${step2}\n"

export assocoutdir=${compoutdir}/gene_lineage_assocations/between_fams_scores

step3="explore network of gene co-evolution (with every single gene lineage represented" 
step3details="(, i.e. with every [collpased clade of] CDS as a single node)"

step4="loading HUGE matrix of co-evolution scores into database"
step4details=", with mention of reconciliation_id to ensure events are not matched across collections"

step5="explore network of gene co-evolution (with network condensed by Orthologous Groups [OGs]" 
step5details=", i.e. with gene lineages belonging to a same OG as a single node)"

if [ ! -z "${somehowusingpostgres}" ] ; then
  echo ${step3}${step3details} "using PostgreSQL db engine"
  ${ptgscripts}/explore_gene_lineage_assoiations.r --dir_match_events ${assocoutdir} --db_type 'postgres'
  checkexec "failed ${step3})" "completed ${step3})\n"
  
  echo ${step4}${step4details}
  ${ptgscripts}/pantagruel_postgres_load_coevolution_scores.py ${sqldb} ${assocoutdir} ${parsedreccolid}
  checkexec "failed ${step4}" "completed ${step4}\n"
  
  echo ${step5}${step5details}
  ${ptgscripts}/explore_orthologous_group_assoiations.r --dir_match_events ${assocoutdir} --db_type 'postgres'
  checkexec "failed ${step5})" "completed ${step5})\n"
  
  python condense_coevolution_network_by_orthologous_groups.py --sqlite_db ${sqldb} 
else
  echo ${step3}${step3details} "using SQLite db engine"
  ${ptgscripts}/explore_gene_lineage_assoiations.r --dir_match_events ${assocoutdir} --db_type 'sqlite'
  checkexec "failed ${step3}" "completed ${step3}\n"
  
  echo """will not load (possibly HUGE) matrix of co-evolution scores into SQLite database, as this is NOT YET IMPLEMENTED to support SQLite db engine. 
  Also, on large data collections this would result into an even larger database file that would be too big to be functional.
  Lack of this prevents the using the orthologous group classification to condense the co-evolution graph by OG (in the current implementation).
  """
  #${ptgscripts}/pantagruel_sqlitedb_load_coevolution_scores.py ${sqldb} ${assocoutdir} ${parsedreccolid}
  exit 0
####
fi