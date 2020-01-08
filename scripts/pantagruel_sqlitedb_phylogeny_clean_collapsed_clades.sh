#!/bin/bash
thisscript=${0}
[ -z "${sqldb}" ] && sqldb="${1}"
[ -z "${collapsecolid}" ] && collapsecolid="${2}"
[ -z "${replacecolid}" ] && replacecolid="${3}"

if [ -z ${replacecolid} ] ; then
 echo "Error: incomplete argument list. Usage:"
 echo "${thisscript} sqlitedb_file collapsed_clade_collection_id replaced_clade_collection_id"
 echo "currently set variable:"
 echo "sqldb=${sqldb} collapsecolid=${collapsecolid} replacecolid=${replacecolid}"
 exit 1
fi

echostep1="delete metadata from collapsed/replaced gene tree collection ${collapsecolid}"
sqlite3 ${sqldb} "DELETE FROM criteria_collapse_gene_tree_clades WHERE criterion_id=${collapsecolid} ;"
sqlite3 ${sqldb} "DELETE FROM criteria_replace_gene_tree_clades WHERE criterion_id=${replacecolid} ;"
checkexec "failed to ${echostep1}" "succeeded to ${echostep1}"


echostep3="delete data of collapsed/replaced gene tree collection ${collapsecolid}"
sqlite3 ${sqldb} """ 
DELETE FROM collapsed_gene_tree_clades WHERE collapse_criterion_id=${collapsecolid};
DELETE FROM replaced_gene_tree_clades WHERE replace_criterion_id=${replacecolid};
"""
checkexec "failed to ${echostep3}" "succeeded to ${echostep3}"

echostep4='delete indexes of collapsed/replaced gene tree collection tables'
sqlite3 ${sqldb} """
DROP INDEX IF EXISTS collapsed_gene_tree_clades_colcritid;
DROP INDEX IF EXISTS collapsed_gene_tree_clades_genefamid;
DROP INDEX IF EXISTS collapsed_gene_tree_clades_cc;
DROP INDEX IF EXISTS collapsed_gene_tree_clades_cdscode ;
DROP INDEX IF EXISTS collapsed_gene_tree_clades_cdscode_colcritid;

DROP INDEX IF EXISTS replaced_gene_tree_clades_replcritid;
DROP INDEX IF EXISTS replaced_gene_tree_clades_genefamid_ccocds;
DROP INDEX IF EXISTS replaced_gene_tree_clades_replab;
DROP INDEX IF EXISTS replaced_gene_tree_clades_replab_replcritid;
"""
checkexec "failed to ${echostep4}" "succeeded to ${echostep4}"