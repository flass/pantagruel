#!/bin/bash
thisscript="${0}"
database="${1}"
dbfile="${2}"
parsedrecs="${3}"
ALEversion="${4}"
ALEalgo="${5}"
ALEsourcenote="${6}"
parsedreccol="${7}"
parsedreccolid="${8}"
parsedreccoldate="${9}"

currsetvar="currently set variable:\ndatabase='${1}' dbfile='${2}' parsedrecs='${3}' ALEversion='${4}' ALEalgo='${5}' ALEsourcenote='${6}' parsedreccol='${7}' parsedreccolid='${8}' parsedreccoldate='${9}'"

if [ -z $parsedreccolid ] ; then
 echo "Error: incomplete argument list. Usage:"
 echo "${thisscript} db_table_dump_folder sqlitedb_file parsed_recs_folder ALE_software_version ALE_bin_name ALE_source_description parsed_recs_collection parsed_recs_collection_id [parsed_recs_collection_date]"
 echo -e ${currsetvar} 
 exit 1
else
 echo -e ${currsetvar} 
fi

cd ${database}

#### NOTE
## here for simplicity only the log variables refering to parsed reconciliations (parsedreccol, parsedreccolid, parsedreccoldate) are recorded in the database
## but not the log variables refering to the actual reconciliations (reccol, reccolid, reccoldate)
####

# load the parsed reconciliation events
python2.7 << EOF 
import sqlite3, os
dbname = '${dbfile}'
dbcon = sqlite3.connect(dbname)
dbcur = dbcon.cursor()
dirgtevt = '${parsedrecs}/gene_tree_lineages'
for nffamgtevt in os.listdir(dirgtevt):
  with open(os.path.join(dirgtevt, nffamgtevt)) as ffamgtevt:
    dbcur.executemany('insert into gene_lineage_events (replacement_label_or_cds_code, event_id, freq, reconciliation_id) values (?,?,?,?);', (line.rstrip('\n').split('\t')+[${parsedreccolid}] for line in ffamgtevt))

dbcon.commit()

dirspet = '${parsedrecs}/ref_species_tree'
with open(os.path.join(dirspet, 'phylogeny_species_tree.tab')) as fspet:
    dbcur.executemany('insert into species_tree (branch_id, parent_branch_id, branch_name, is_tip) values (?,?,?,?);', (line.rstrip('\n').split('\t') for line in fspet))

with open(os.path.join(dirspet, 'phylogeny_species_tree_events.tab')) as fspetevt:
    dbcur.executemany('insert into species_tree_events (event_id, event_type, don_branch_id, rec_branch_id) values (?,?,?,?);', (line.rstrip('\n').split('\t') for line in fspetevt))

dbcon.commit()
EOF
# consider using the following for speeding up:
#~ df = pandas.read_csv(csvfile)
#~ df.to_sql(table_name, conn, if_exists='append', index=False)

sqlite3 ${dbfile} << EOF 
BEGIN;
INSERT INTO reconciliation_collections (reconciliation_id, reconciliation_name, software, version, algorithm, reconciliation_date, notes)
 VALUES (${parsedreccolid}, '${parsedreccol}', 'ALE', '${ALEversion}', '${ALEalgo}', '${parsedreccoldate}', '${ALEsourcenote}') ;

CREATE INDEX gene_lineage_events_recid ON gene_lineage_events (reconciliation_id);
CREATE INDEX gene_lineage_events_rlocds ON gene_lineage_events (replacement_label_or_cds_code);
CREATE INDEX gene_lineage_events_evtid ON gene_lineage_events (event_id);
CREATE INDEX gene_lineage_events_freq ON gene_lineage_events (freq);
CREATE INDEX gene_lineage_events_rlocds_evtid ON gene_lineage_events (replacement_label_or_cds_code, event_id);
CREATE UNIQUE INDEX gene_lineage_events_recid_rlocds_evtid ON gene_lineage_events (reconciliation_id, replacement_label_or_cds_code, event_id);

INSERT INTO replacement_label_or_cds_code2gene_families (replacement_label_or_cds_code, gene_family_id) 
SELECT replacement_label_or_cds_code, gene_family_id FROM (
 SELECT cds_code as replacement_label_or_cds_code, gene_family_id FROM coding_sequences
UNION
 SELECT replacement_label as replacement_label_or_cds_code, gene_family_id FROM replaced_gene_tree_clades
) q1 
INNER JOIN (SELECT DISTINCT replacement_label_or_cds_code FROM gene_lineage_events) q2 USING (replacement_label_or_cds_code);

CREATE UNIQUE INDEX rlocds2genefam_rlocds ON replacement_label_or_cds_code2gene_families (replacement_label_or_cds_code);
CREATE INDEX rlocds2genefam_genefam ON replacement_label_or_cds_code2gene_families (gene_family_id);

INSERT INTO gene_tree_label2cds_code (replacement_label_or_cds_code, cds_code) 
SELECT replacement_label_or_cds_code, cds_code FROM (
  SELECT cds_code as replacement_label_or_cds_code, cds_code 
   FROM coding_sequences
 UNION 
  SELECT rgtc1.replacement_label as replacement_label_or_cds_code, rgtc1.col_clade_or_cds_code as cds_code
   FROM replaced_gene_tree_clades AS rgtc1
   WHERE rgtc1.col_clade_or_cds_code NOT LIKE 'clade%'
 UNION
  SELECT rgtc2.replacement_label as replacement_label_or_cds_code, cgtc.cds_code
   FROM replaced_gene_tree_clades AS rgtc2
   INNER JOIN collapsed_gene_tree_clades AS cgtc ON rgtc2.col_clade_or_cds_code=cgtc.col_clade AND rgtc2.gene_family_id=cgtc.gene_family_id
) q1 
INNER JOIN (SELECT DISTINCT replacement_label_or_cds_code FROM gene_lineage_events) q2 USING (replacement_label_or_cds_code)
;

CREATE INDEX gtlab2cds_rlocds ON gene_tree_label2cds_code (replacement_label_or_cds_code);
CREATE UNIQUE INDEX gtlab2cds_cdscode ON gene_tree_label2cds_code (cds_code);

COMMIT;
VACUUM;
EOF
