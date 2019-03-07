#!/bin/bash
thisscript=${0}
database=${1}
dbfile=${2}
colalinexuscodedir=${3}
coltreechains=${4}
collapsecond=${5}
colmethod=${6}
collapsecriteriondef=${7}
collapsecolid=${8}
replacecolid=${9}
collapsecoldate=${10}
replacecoldate=${11}

if [ -z $database ] ; then
 echo "Error: incomplete argument list. Usage:"
 echo "${thisscript} db_table_dump_folder sqlitedb_file collapsed_alignment_folder collapsed_treechain_folder collapse_criterion replacement_criterion [collapsed_clade_collection_id replaced_clade_collection_id collapsed_clade_collection_date replaced_clade_collection_date]"
 echo "currently set variable:"
 echo "database=$1 dbfile=$2 colalinexuscodedir=$3 coltreechains=$4 collapsecond=$5 colmethod=$6 collapsecriteriondef=$7 collapsecolid=$8 replacecolid=$9 collapsecoldate=$10 replacecoldate=$11"
 exit 1
fi
if [ -z ${collapsecolid} ] ; then
  collapsecolid=0
fi
if [ -z ${replacecolid} ] ; then
  replacecolid=0
fi
cd ${database}

# load data on collapsed gene tree clades and their replacement in reconciled gene trees
sqlite3 ${dbfile} << EOF 
INSERT INTO criteria_collapse_gene_tree_clades (criterion_id, criterion_name, criterion_definition, collapsed_clade_collection_creation)
 VALUES (${collapsecolid}, '${collapsecond}', '$(echo ${collapsecriteriondef} | sed -e "s/'/''/g")', '${collapsecoldate}') ;
INSERT INTO criteria_replace_gene_tree_clades (criterion_id, criterion_name, criterion_definition, replaced_clade_collection_creation)
 VALUES (${replacecolid}, '${colmethod}', '--method=${colmethod}', '${replacecoldate}') ;
EOF

# create big table files for collapsed gene tree clades and replaced clades in reconciled gene trees
python ${ptgscripts}/generate_collapsed_replaced_clades_tables.py \
 ${colalinexuscodedir}/${collapsecond}/mbconstraints ${coltreechains}/${collapsecond}/${colmethod} \
 ${collapsecolid} ${replacecolid} \
 ${database}/phylogeny_collapsed_gene_tree_clades.tab ${database}/phylogeny_replaced_gene_tree_clades.tab
# load the in db
sqlite3 ${dbfile} << EOF 
.mode tabs
.import '${database}/phylogeny_collapsed_gene_tree_clades.tab' collapsed_gene_tree_clades
.import '${database}/phylogeny_replaced_gene_tree_clades.tab' replaced_gene_tree_clades
EOF
