#!/bin/bash
thisscript=$0
database=$1
dbfile=$2

if [ -z $cdsorfanclust ] ; then
 echo "Error: incomplete argument list. Usage:"
 echo "${thisscript} db_table_dump_folder sqlitedb_file"
 echo "currently set variable:"
 echo "database=$1 dbfile=$2"
 exit 1
fi


initiatescript=${thisscript%.*}_initiate.sql
populatescript=${thisscript%.*}_populate.py

cd ${database}

### if no set db name env variable, INTERACTIVE prompt asks for database name and password
while [ -z $dbfile ] ; do
  read -p 'Set PostgreSQL database name: ' dbfile ; echo ''
done

# load database (phylogeny schema)
sqlite3 ${dbfile} < $initiatescript

sqlite3 ${dbfile} < $ptgscripts/pantagruel_postgres_phylogeny_schema.sql
# load data on collapsed gene tree clades and their replacement in reconciled gene trees
sqlite3 ${dbfile} << EOF 
INSERT INTO criteria_collapse_gene_tree_clades (criterion_id, criterion_name, criterion_definition, collapsed_clade_collection_creation)
 VALUES (1, 'bs_stem70_withinmedian35', E'--clade_stem_conds="[(\'bs\', \'>=\', 70)]" --within_clade_conds="[(\'median\', \'bs\', \'<=\', 35, -1), (\'max\', \'bs\', \'<\', 70, -1)]"', '2017-11-15') ;
INSERT INTO criteria_replace_gene_tree_clades (criterion_id, criterion_name, criterion_definition, replaced_clade_collection_creation)
 VALUES (1, 'replaceCCinGasinS-collapsePOPinSnotinG', '', '2018-02-20') ;
EOF
# create big table files for collapsed gene tree clades and replaced clades in reconciled gene trees
python $ptgscripts/generate_collapsed_replaced_clades_tables.py \
 ${colalinexuscodedir}/${collapsecond}/mbconstraints ${coltreechains}/${collapsecond}/${colmethod} \
 ${database}/phylogeny_collapsed_gene_tree_clades.tab ${database}/phylogeny_replaced_gene_tree_clades.tab 
# load the in db
sqlite3 ${dbfile} << EOF 
.mode tabs
.import '${database}/phylogeny_collapsed_gene_tree_clades.tab collapsed_gene_tree_clades (gene_family_id, col_clade, cds_code)'
UPDATE collapsed_gene_tree_clades SET collapse_criterion_id=1 WHERE collapse_criterion_id IS NULL;
.import '${database}/phylogeny_replaced_gene_tree_clades.tab replaced_gene_tree_clades (gene_family_id, col_clade_or_cds_code, replacement_label)'
UPDATE replaced_gene_tree_clades SET replace_criterion_id=1 WHERE replace_criterion_id IS NULL;
EOF
# load the parsed reconciliation events
python << EOF 
import sqlite3, os
dbname = '${dbfile}'
dbcon = sqlite3.connect(dbname=dbname)
dbcur = dbcon.cursor()
dirgtevt = '${compoutdir}/ale_collapsed_undat/gene_tree_lineages'
for nffamgtevt in os.listdir(dirgtevt):
  with open(os.path.join(dirgtevt, nffamgtevt)) as ffamgtevt:
    dbcur.executemany('insert into gene_lineage_events (replacement_label_or_cds_code, event_id, freq) vaues (?,?,?);', (line.rstrip('\n').split('\t') for line in ffamgtevt))

dbcon.commit()
EOF
# consider using the following for speeding up:
#~ df = pandas.read_csv(csvfile)
#~ df.to_sql(table_name, conn, if_exists='append', index=False)

sqlite3 ${dbfile} << EOF 

INSERT INTO reconciliation_collections (reconciliation_id, reconciliation_name, software, version, algorithm, reconciliation_date, notes)
 VALUES (1, 'ale_collapsed_undat', 'ALE', 'v0.4', 'ALEml_undated', '2018-02-28', E'program compiled from source code from branch \'master\' of https://github.com/ssolo/ALE commit 63f0a3c964074a15f61fd45156ab9e10b5dd45ef') ;
UPDATE gene_lineage_events SET reconciliation_id=1 WHERE reconciliation_id IS NULL;

CREATE INDEX ON gene_lineage_events (reconciliation_id);
CREATE INDEX ON gene_lineage_events (replacement_label_or_cds_code);
CREATE INDEX ON gene_lineage_events (event_id);
CREATE INDEX ON gene_lineage_events (freq);
ALTER TABLE gene_lineage_events ADD PRIMARY KEY (event_id, cds_code, reconciliation_id);

CREATE INDEX ON collapsed_gene_tree_clades (gene_family_id);
CREATE INDEX ON collapsed_gene_tree_clades (cds_code);
CREATE INDEX ON collapsed_gene_tree_clades (gene_family_id, col_clade);

CREATE INDEX ON replaced_gene_tree_clades (gene_family_id, col_clade_or_cds_code);
CREATE INDEX ON replaced_gene_tree_clades (replacement_label);

CREATE TABLE replacement_label_or_cds_code2gene_families AS 
SELECT cds_code as replacement_label_or_cds_code, gene_family_id FROM coding_sequences
UNION
SELECT replacement_label as replacement_label_or_cds_code, gene_family_id FROM replaced_gene_tree_clades;
DELETE FROM replacement_label_or_cds_code2gene_families
 WHERE replacement_label_or_cds_code NOT IN (SELECT DISTINCT replacement_label_or_cds_code FROM gene_lineage_events);
-- ALTER TABLE replacement_label_or_cds_code2gene_families ADD COLUMN rlocds_id SERIAL PRIMARY KEY;
-- use OIDs instead

CREATE TABLE gene_tree_label2cds_code (replacement_label_or_cds_code, cds_code) AS
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
;
DELETE FROM gene_tree_label2cds_code
 WHERE replacement_label_or_cds_code NOT IN (SELECT DISTINCT replacement_label_or_cds_code FROM gene_lineage_events);

CREATE INDEX ON replacement_label_or_cds_code2gene_families (replacement_label_or_cds_code);
CREATE INDEX ON replacement_label_or_cds_code2gene_families (gene_family_id);
CREATE INDEX ON gene_tree_label2cds_code (replacement_label_or_cds_code);
CREATE INDEX ON gene_tree_label2cds_code (cds_code);

VACUUM ;
EOF
