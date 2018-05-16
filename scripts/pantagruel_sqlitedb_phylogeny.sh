#!/bin/bash
thisscript=$0
database=$1
dbname=$2

if [ -z $cdsorfanclust ] ; then
 echo "Error: incomplete argument list. Usage:"
 echo "${thisscript} database dbname"
 echo "currently set variable:"
 echo "database=$1 dbname=$2"
 exit 1
fi


initiatescript=${thisscript%.*}_initiate.sql
populatescript=${thisscript%.*}_populate.py

cd ${database}

### if no set db name env variable, INTERACTIVE prompt asks for database name and password
while [ -z $dbname ] ; do
  read -p 'Set PostgreSQL database name: ' dbname ; echo ''
done

# load database (phylogeny schema)
sqlite3 ${dbname} < $initiatescript

sqlite3 ${dbname} < $ptgscripts/pantagruel_postgres_phylogeny_schema.sql
# load data on collapsed gene tree clades and their replacement in reconciled gene trees
sqlite3 ${dbname} << EOF 
INSERT INTO phylogeny.criteria_collapse_gene_tree_clades (criterion_id, criterion_name, criterion_definition, collapsed_clade_collection_creation)
 VALUES (1, 'bs_stem70_withinmedian35', E'--clade_stem_conds="[(\'bs\', \'>=\', 70)]" --within_clade_conds="[(\'median\', \'bs\', \'<=\', 35, -1), (\'max\', \'bs\', \'<\', 70, -1)]"', '2017-11-15') ;
INSERT INTO phylogeny.criteria_replace_gene_tree_clades (criterion_id, criterion_name, criterion_definition, replaced_clade_collection_creation)
 VALUES (1, 'replaceCCinGasinS-collapsePOPinSnotinG', '', '2018-02-20') ;
EOF
# create big table files for collapsed gene tree clades and replaced clades in reconciled gene trees
python $ptgscripts/generate_collapsed_replaced_clades_tables.py \
 ${colalinexuscodedir}/${collapsecond}/mbconstraints ${coltreechains}/${collapsecond}/${colmethod} \
 ${database}/phylogeny_collapsed_gene_tree_clades.tab ${database}/phylogeny_replaced_gene_tree_clades.tab 
# load the in db
sqlite3 ${dbname} << EOF 
\copy phylogeny.collapsed_gene_tree_clades (gene_family_id, col_clade, cds_code) from '${database}/phylogeny_collapsed_gene_tree_clades.tab'
UPDATE phylogeny.collapsed_gene_tree_clades SET collapse_criterion_id=1 WHERE collapse_criterion_id IS NULL;
\copy phylogeny.replaced_gene_tree_clades (gene_family_id, col_clade_or_cds_code, replacement_label) from '${database}/phylogeny_replaced_gene_tree_clades.tab'
UPDATE phylogeny.replaced_gene_tree_clades SET replace_criterion_id=1 WHERE replace_criterion_id IS NULL;
EOF
# load the parsed reconciliation events
python << EOF 
import sqlite3, os
dbname = '${sqldbname}'
dbcon = sqlite3.connect(dbname=dbname)
dbcur = dbcon.cursor()
dirgtevt = '${compoutdir}/ale_collapsed_undat/gene_tree_lineages'
for nffamgtevt in os.listdir(dirgtevt):
  with open(os.path.join(dirgtevt, nffamgtevt)) as ffamgtevt:
    dbcur.executemany('insert into phylogeny.gene_lineage_events (cds_code, event_id, freq) vaues (?,?,?);', (line.rstrip('\n').split('\t') for line in ffamgtevt))

dbcon.commit()
EOF
consider using the following for speeding up:
#~ df = pandas.read_csv(csvfile)
#~ df.to_sql(table_name, conn, if_exists='append', index=False)


