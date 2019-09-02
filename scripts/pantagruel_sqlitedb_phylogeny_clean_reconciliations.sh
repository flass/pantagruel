#!/bin/bash
thisscript="${0}"
database="${1}"
dbfile="${2}"
parsedreccolid="${3}"

currsetvar="currently set variable:\ndatabase='${1}' dbfile='${2}' parsedreccolid='${3}'"

if [ -z $parsedreccolid ] ; then
 echo "Error: incomplete argument list. Usage:"
 echo "${thisscript} db_table_dump_folder sqlitedb_file parsed_recs_collection_id"
 echo -e ${currsetvar} 
 exit 1
else
 echo -e ${currsetvar} 
fi

cd ${database}


sqlite3 ${dbfile} << EOF 
BEGIN;

DELETE FROM species_tree;
DELETE FROM species_tree_events;
DELETE FROM gene_lineage_events WHERE reconciliation_id=${parsedreccolid};

DELETE FROM reconciliation_collections WHERE reconciliation_id=${parsedreccolid};


DROP INDEX IF EXISTS gene_lineage_events_recid;
DROP INDEX IF EXISTS gene_lineage_events_rlocds;
DROP INDEX IF EXISTS gene_lineage_events_evtid;
DROP INDEX IF EXISTS gene_lineage_events_freq;
DROP INDEX IF EXISTS gene_lineage_events_rlocds_evtid;
DROP INDEX IF EXISTS ;

DROP INDEX IF EXISTS collapsed_gene_tree_clades_colcritid;
DROP INDEX IF EXISTS collapsed_gene_tree_clades_genefamid;
DROP INDEX IF EXISTS collapsed_gene_tree_clades_cc;
DROP INDEX IF EXISTS collapsed_gene_tree_clades_cdscode;
DROP INDEX IF EXISTS collapsed_gene_tree_clades_cdscode_colcritid;

DROP INDEX IF EXISTS replaced_gene_tree_clades_replcritid;
DROP INDEX IF EXISTS replaced_gene_tree_clades_genefamid_ccocds;
DROP INDEX IF EXISTS replaced_gene_tree_clades_replab;
DROP INDEX IF EXISTS replaced_gene_tree_clades_replab_replcritid;

DELETE FROM replacement_label_or_cds_code2gene_families ;

DROP UNIQUE INDEX rlocds2genefam_rlocds;
DROP INDEX IF EXISTS rlocds2genefam_genefam;

DELETE FROM gene_tree_label2cds_code;

DROP INDEX IF EXISTS gtlab2cds_rlocds;
DROP UNIQUE INDEX gtlab2cds_cdscode;

COMMIT;
VACUUM;
EOF
