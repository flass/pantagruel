#!/bin/bash
thisscript=${0}
[ -z "${database}" ] && database="${1}"
[ -z "${sqldb}" ] && sqldb="${2}"
[ -z "${colalinexuscodedir}" ] && colalinexuscodedir="${3}"
[ -z "${coltreechains}" ] && coltreechains="${4}"
[ -z "${collapsecond}" ] && collapsecond="${5}"
[ -z "${replmethod}" ] && replmethod="${6}"
[ -z "${collapsecriteriondef}" ] && collapsecriteriondef="${7}"
[ -z "${collapsecolid}" ] && collapsecolid="${8}"
[ -z "${replacecolid}" ] && replacecolid="${9}"
[ -z "${collapsecoldate}" ] && collapsecoldate="${10}"
[ -z "${replacecoldate}" ] && replacecoldate="${11}"

if [ -z $database ] ; then
 echo "Error: incomplete argument list. Usage:"
 echo "${thisscript} db_table_dump_folder sqlitedb_file collapsed_alignment_folder collapsed_treechain_folder collapse_criterion replacement_criterion [collapsed_clade_collection_id replaced_clade_collection_id collapsed_clade_collection_date replaced_clade_collection_date]"
 echo "currently set variable:"
 echo "database=${1} sqldb=${2} colalinexuscodedir=${3} coltreechains=${4} collapsecond=${5} replmethod=${6} collapsecriteriondef=${7} collapsecolid=${8} replacecolid=${9} collapsecoldate=${10} replacecoldate=${11}"
 exit 1
fi
if [ -z ${collapsecolid} ] ; then
  collapsecolid=0
fi
if [ -z ${replacecolid} ] ; then
  replacecolid=0
fi
cd ${database}

dblquotecollapsecriteriondef=$(echo ${collapsecriteriondef} | sed -e "s/'/''/g")
echostep1='load metadata on (collapsed) gene tree collections'
sqlite3 ${sqldb} << EOF 
INSERT INTO criteria_collapse_gene_tree_clades (criterion_id, criterion_name, criterion_definition, collapsed_clade_collection_creation)
 VALUES (${collapsecolid}, '${collapsecond}', '${dblquotecollapsecriteriondef}', '${collapsecoldate}') ;
INSERT INTO criteria_replace_gene_tree_clades (criterion_id, criterion_name, criterion_definition, replaced_clade_collection_creation)
 VALUES (${replacecolid}, '${replmethod}', '--method=${replmethod}', '${replacecoldate}') ;
EOF
checkexec "failed to ${echostep1}" "succeeded to ${echostep1}"

echostep2='create big table files for collapsed gene tree clades and their replacement clades in tree chains to be reconciled'
python2.7 ${ptgscripts}/generate_collapsed_replaced_clades_tables.py \
 ${colalinexuscodedir}/${collapsecond}/mbconstraints ${coltreechains}/${collapsecond}/${replmethod} \
 ${collapsecolid} ${replacecolid} \
 ${database}/phylogeny_collapsed_gene_tree_clades.tab ${database}/phylogeny_replaced_gene_tree_clades.tab
checkexec "failed to ${echostep2}" "succeeded to ${echostep2}"

echostep3='load big table files of collapsed/replacement clade data into the sqlite db'
sqlite3 ${sqldb} << EOF
.mode tabs
.import '${database}/phylogeny_collapsed_gene_tree_clades.tab' collapsed_gene_tree_clades
.import '${database}/phylogeny_replaced_gene_tree_clades.tab' replaced_gene_tree_clades
EOF
checkexec "failed to ${echostep3}" "succeeded to ${echostep3}"


echostep4='index tables of collapsed/replacement clade data into the sqlite db'
sqlite3 ${sqldb} """
BEGIN;
CREATE INDEX collapsed_gene_tree_clades_colcritid ON collapsed_gene_tree_clades (collapse_criterion_id);
CREATE INDEX collapsed_gene_tree_clades_genefamid ON collapsed_gene_tree_clades (gene_family_id);
CREATE INDEX collapsed_gene_tree_clades_cc ON collapsed_gene_tree_clades (gene_family_id, col_clade);
CREATE INDEX collapsed_gene_tree_clades_cdscode ON collapsed_gene_tree_clades (cds_code);
CREATE UNIQUE INDEX collapsed_gene_tree_clades_cdscode_colcritid ON collapsed_gene_tree_clades (cds_code, collapse_criterion_id);

CREATE INDEX replaced_gene_tree_clades_replcritid ON replaced_gene_tree_clades (replace_criterion_id);
CREATE INDEX replaced_gene_tree_clades_genefamid_ccocds ON replaced_gene_tree_clades (gene_family_id, col_clade_or_cds_code);
CREATE INDEX replaced_gene_tree_clades_replab ON replaced_gene_tree_clades (replacement_label);
CREATE UNIQUE INDEX replaced_gene_tree_clades_replab_replcritid ON replaced_gene_tree_clades (replacement_label, replace_criterion_id);
COMMIT;
ANALYZE;
"""
checkexec "failed to ${echostep4}" "succeeded to ${echostep4}"
