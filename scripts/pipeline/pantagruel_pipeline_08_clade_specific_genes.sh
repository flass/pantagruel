#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

if [ -z "$2" ] ; echo "missing mandatory parameters." ; echo "Usage: $0 ptg_db_name ptg_root_folder" ; exit 1 ; fi
export ptgdbname="$1"  # database anme (will notably be the name of the top folder)
export ptgroot="$2"    # source folder where to create the database
envsourcescript=${ptgroot}/environ_pantagruel_${ptgdbname}.sh
source $envsourcescript

###############################################
## 08. orthologous and clade-specific gene sets
###############################################

export orthogenes=${ptgdb}/08.orthologs
mkdir -p ${orthogenes}

cd ${ptgrepo} ; export ptgversion=$(git log | grep commit | cut -d' ' -f2) ; cd -

# classify genes into orthologous groups for each gene of the reconciled gene tree sample
# do not report detailed results but, using gptgh analysis, combine the sample-wide classification into one classification for the gene family
# run in parallel

## for the moment only coded for dated ALE model (ALEml reconciliations)
## and with no colpased gene tree clades (need to discard events below replacement clade subtree roots [as in parse_collapsedALE_scenarios.py]
## and to transpose orthologus group classification of collapsed clade to member genes)

if [ -z ${getOrpthologuesOptions} ] ; then
  getOrpthologuesOptions=" --ale.model='dated' --methods='mixed' --max.frac.extra.spe=0.5 --majrule.combine=0.5 --colour.combined.tree"
fi
if [ -z ${orthoColId} ] ; then
  orthocolid=1
fi
# generate Ortholog Collection
orthocol=ortholog_collection_${orthoColId}
mkdir -p ${orthogenes}/${orthocol}
${ptgscripts}/get_orthologues_from_ALE_recs.py -i ${outrecdir} -o ${orthogenes}/${orthocol} ${getOrpthologuesOptions} -T 8 &> $ptglogs/get_orthologues_from_ALE_recs_${orthocol}.log

# import ortholog classification into database
sqlite3 ${sqldb} """INSERT INTO ortholog_collections (ortholog_col_id, ortholog_col_name, reconciliation_id, software, version, algorithm, ortholog_col_date, notes) VALUES 
(${orthocolid}, '${orthocol}', ${parsedreccolid}, 'pantagruel/scripts/get_orthologues_from_ALE_recs.py', '${ptgversion:0:7}', 'getOrthologues(method=''mixed'')', '$(date +%Y-%m-%d)', 
'source from https://github.com/flass/pantagruel/commits/${ptgversion}, call: ''scripts/get_orthologues_from_ALE_recs.py ${getOrpthologuesOptions}''')
;
"""
python ${ptgscripts}/pantagruel_sqlitedb_load_orthologous_groups.py ${sqldb} ${orthogenes}/ortholog_collection_1 "mixed" "majrule_combined_0.500000" ${orthocolid}


# generate abs/pres matrix
orthocol=ortholog_collection_${orthocolid}
echo $orthocol
orthomatrad=${orthogenes}/${orthocol}/mixed_majrule_combined_0.5.orthologs
python ${ptgscripts}/get_ortholog_presenceabsence_matrix_from_sqlitedb.py ${sqldb} ${orthomatrad} ${coregenome}/${focus}/${focus}_genome_codes ${orthocolid}
${ptgscripts}/get_clade_specific_genes.r ${orthomatrad}_genome_counts.no-singletons.mat ${sqldb} ${orthocolid} ${coregenome}/${focus}/${focus} ${orthomatrad}

# list clade-specific orthologs
export orthomatrad=${orthogenes}/${orthocol}/mixed_majrule_combined_0.5.orthologs
