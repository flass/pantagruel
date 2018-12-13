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
export ptgdb=${ptgroot}/${ptgdbname}
envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
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
  getOrpthologuesOptions=" --ale.model=${rectype} --methods='mixed' --max.frac.extra.spe=0.5 --majrule.combine=0.5 --colour.combined.tree --use.unreconciled.gene.trees ${mbgenetrees} --unreconciled.format 'nexus' --unreconciled.ext '.con.tre'"
fi
if [ -z ${orthoColId} ] ; then
  orthocolid=1
fi
# derived parameters
if [ ${ALEalgo} == 'ALEml_undated' ] ; then
  export rectype='undat'
else
  export rectype='dated'
fi
# generate Ortholog Collection
orthocol=ortholog_collection_${orthoColId}
mkdir -p ${orthogenes}/${orthocol}
${ptgscripts}/get_orthologues_from_ALE_recs.py -i ${outrecdir} -o ${orthogenes}/${orthocol} ${getOrpthologuesOptions} &> $ptglogs/get_orthologues_from_ALE_recs_${orthocol}.log

# import ortholog classification into database
sqlite3 ${sqldb} """INSERT INTO ortholog_collections (ortholog_col_id, ortholog_col_name, reconciliation_id, software, version, algorithm, ortholog_col_date, notes) VALUES 
(${orthocolid}, '${orthocol}', ${parsedreccolid}, 'pantagruel/scripts/get_orthologues_from_ALE_recs.py', '${ptgversion:0:7}', 'getOrthologues(method=''mixed'')', '$(date +%Y-%m-%d)', 
'source from https://github.com/flass/pantagruel/commits/${ptgversion}, call: ''scripts/get_orthologues_from_ALE_recs.py ${getOrpthologuesOptions}''')
;
"""
python ${ptgscripts}/pantagruel_sqlitedb_load_orthologous_groups.py ${sqldb} ${orthogenes}/${orthocol} "mixed" "majrule_combined_0.500000" ${orthocolid}
python ${ptgscripts}/pantagruel_sqlitedb_load_orthologous_groups.py ${sqldb} ${orthogenes}/${orthocol} "unreconciled" "unicopy" ${orthocolid} 1

## extract sister clade pairs from the reference species tree for later clade-specific gene search
export cladedefs="${speciestree}_clade_defs"
python << EOF
import tree2
reftree = tree2.Node(file='${speciestree}')
nfout = '${cladedefs}'
fout = open(nfout, 'w')
fout.write('\t'.join(['', 'clade', 'sisterclade', 'name', 'maxabsin', 'maxpresout'])+'\n')
k = 0
for node in reftree:
  if len(node.children) != 2: continue
  for child in node.children:
    if child.nb_leaves() <= 1: continue
    focchildlab = child.label()
    if not focchildlab:
      focchildlab = "clade%d"%k
      k += 1
    focchildleaflabset = ','.join(sorted(child.get_leaf_labels()))
    sischildleaflabset = ','.join(sorted(child.go_brother().get_leaf_labels()))
    fout.write('\t'.join([focchildlab, '""', focchildleaflabset, sischildleaflabset, '0', '0'])+'\n')

fout.close()
EOF
# the resulting fle ${cladedefs} can have information filled/changed in the following collumns:
# - 'name':       a taxon or practical name for the clade
# - 'maxabsin':   the maximum count of genomes in the focal clade missing the gene (default 0)
# - 'maxpresout': the maximum count of genomes in the sister clade featuring the gene (default 0)
# additionally, extra rows can be added with custom genome group definitions (they don't need to be clades)

# generate abs/pres matrix
orthocol=ortholog_collection_${orthocolid}
echo ${orthocol}
orthomatrad=${orthogenes}/${orthocol}/mixed_majrule_combined_0.5.orthologs
python ${ptgscripts}/get_ortholog_presenceabsence_matrix_from_sqlitedb.py ${sqldb} ${orthomatrad} ${orthocolid}

# create clsutering based on the abs/pres matrix (using Jaccard Distance)
${dbscripts}/pangenome_hclust.r ${orthomatrad} 1000 & 

# list clade-specific orthologs
${ptgscripts}/get_clade_specific_genes.r ${orthomatrad}_genome_counts.no-singletons.mat ${sqldb} ${orthocolid} ${speciestree} ${orthomatrad}

## test GO term enrichment in gene sets

## generate background term distribution for clades i.e. list all genes in their pangenome and associated terms
export goterms=${funcannot}/GeneOntology
mkdir -p ${goterms}
# for the whole dataset
sqlite3 -cmd ".mode tab" ${sqldb} "select distinct locus_tag, go_id from coding_sequences 
left join functional_annotations using (nr_protein_id) 
left join interpro2GO using (interpro_id) ;" > ${goterms}/${ngenomes}-genomes_pangenome_terms.tab
# for all clades of the species tree
export claderefgodir=${goterms}/clade_go_term_reference_sets
mkdir -p ${claderefgodir}/
tail -n +2 ${cladedefs} | while read cla cladedef siscladedef ; do
  claspeset="'$(echo $cladedef | sed -e "s/,/','/g")'"
  echo $cla $claspeset
  cladest=${claderefgodir}/${cla}_pangenome_terms.tab
  q="select distinct locus_tag, go_id from coding_sequences 
  inner join replicons using (genomic_accession)
  inner join assemblies using (assembly_id)
  left join functional_annotations using (nr_protein_id) 
  left join interpro2GO using (interpro_id) 
  where code in (${claspeset})"
  sqlite3 -cmd ".mode tab" ${sqldb} "${q};" > ${cladest}
  sqlite3 -cmd ".mode tab" ${sqldb} "${q} and go_id not null;" > ${cladest}_nonull
  ls -lh ${cladest}
done

export dirgotablescladespe=${orthomatrad}_specific_genes.tables_byclade_goterms_pathways
export dirgoenrichcladespe=${goterms}/clade_go_term_enriched_cladespecific_vs_coregenome
mkdir -p ${dirgoenrichcladespe}/
# compare each clade-specific core genome (single repr sequence per OG) to its respective core genome (single repr sequence per OG)
tail -n +2 ${cladedefs} | while read cla name cladedef siscladedef maxabs maxpres ; do
  echo $cla
  cladspego=${dirgotablescladespe}/mixed_majrule_combined_0.5.orthologs_specific_genes_${cla}_reprseq_goterms.tab
  cut -f5,6 ${cladspego} | grep -v "NA$" > ${cladspego}_nonull
  ${dbscripts}/clade_specific_genes_GOterm_enrichment_test.r \
  --study_annots ${cladspego}_nonull  \
  --population_annots ${claderefgodir}/${cla}_coregenome_terms.tab_nonull \
  --out ${dirgoenrichcladespe}/${cla}_go_term_enriched_cladespecific_vs_coregenome.tab \
  --algo "weight01" --stat "Fisher" &> ${raplogs}/clade_specific_vs_coregenome_${cla}_GOterm_enrichment_test.log
  ls -lh ${dirgoenrichcladespe}/*_${cla}_* ; echo ""
done &> ${raplogs}/cladespecific_vs_coregenome_genes_GOterm_enrichment_test.log

export dirgotablescladespe=${orthomatrad}_specific_genes.tables_byclade_goterms_pathways
export dirgoenrichcladespe=${goterms}/clade_go_term_enriched_cladespecific_vs_pangenome
mkdir -p ${dirgoenrichcladespe}/
# compare each clade-specific core genome (all sequences) to its respective pangenome (all sequences)
tail -n +2 ${cladedefs} | while read cla name cladedef siscladedef maxabs maxpres ; do
  echo $cla
  cladspego=${dirgotablescladespe}/mixed_majrule_combined_0.5.orthologs_specific_genes_${cla}_allseq_goterms.tab
  cut -f5,6 ${cladspego} | grep -v "NA$" > ${cladspego}_nonull
  ${dbscripts}/clade_specific_genes_GOterm_enrichment_test.r \
  --study_annots ${cladspego}_nonull  \
  --population_annots ${claderefgodir}/${cla}_pangenome_terms.tab_nonull \
  --out ${dirgoenrichcladespe}/${cla}_go_term_enriched_cladespecific_vs_pangenome.tab \
  --algo "weight01" --stat "Fisher" &> ${raplogs}/clade_specific_vs_coregenome_${cla}_GOterm_enrichment_test.log
  ls -lh ${dirgoenrichcladespe}/*${cla}* ; echo ""
done &> ${raplogs}/cladespecific_vs_pangenome_genes_GOterm_enrichment_test.log

# concatenate summary reports
tail -n +2 ${cladedefs} | while read cla name cladedef siscladedef maxabs maxpres ; do
  cla=clade${n}
  echo "# ${cla} ${name}"
  echo "# cladespe vs. core"
  cat ${goterms}/41NeoPseudo_clade_go_term_enriched_cladespecific_vs_coregenome/41NeoPseudo_${cla}_go_term_enriched_cladespecific_vs_coregenome.tab
  echo "# cladespe vs. pan"
  cat ${goterms}/41NeoPseudo_clade_go_term_enriched_cladespecific_vs_pangenome/41NeoPseudo_${cla}_go_term_enriched_cladespecific_vs_pangenome.tab
  echo "# - - - "
done > ${goterms}/41NeoPseudo_clade_go_term_enriched_cladespecific_summary.tab
