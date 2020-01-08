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

echostep1="delete metadata from (collapsed) gene tree collection ${collapsecolid}"
sqlite3 ${sqldb} "DELETE FROM criteria_collapse_gene_tree_clades WHERE criterion_id=${collapsecolid} ;"
checkexec "failed to ${echostep1}" "succeeded to ${echostep1}"


echostep3="delete data of (collapsed) gene tree collection ${collapsecolid}"
sqlite3 ${sqldb} """ 
DELETE FROM collapsed_gene_tree_clades WHERE collapse_criterion_id=${collapsecolid};
DELETE FROM replaced_gene_tree_clades WHERE replace_criterion_id=${replacecolid};
"""
checkexec "failed to ${echostep3}" "succeeded to ${echostep3}"
