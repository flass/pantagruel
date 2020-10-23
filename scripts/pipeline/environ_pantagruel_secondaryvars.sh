#!/usr/bin/env bash

## This source file should be called by a Pantagruel database configuration file 
## to generate all database-related secondary environment variables

# head folders
export ptgdb=${ptgroot}/${ptgdbname}
export ptglogs=${ptgdb}/logs
export ptgtmp=${ptgdb}/tmp
export indata=${ptgdb}/00.input_data
export seqdb=${ptgdb}/01.seqdb
export protali=${ptgdb}/02.gene_alignments
export database=${ptgdb}/03.database
export funcannot=${ptgdb}/04.functional
export coregenome=${ptgdb}/05.core_genome
export genetrees=${ptgdb}/06.gene_trees
export alerec=${ptgdb}/07.reconciliations
export orthogenes=${ptgdb}/08.orthologs
export comparerecs=${ptgdb}/09.compare_scenarios

# sub folders
export contigs=${customassemb}/contigs
export custannot=${customassemb}/annotation
export annot=${indata}/annotation
export prokkaref=${indata}/reference_assemblies4annotation
export gblikeass=${indata}/genbank-format_assemblies
export genomeinfo=${indata}/genome_infos
export assemblies=${indata}/assemblies
export gp2ass=${indata}/genomesource_assemblyid_assemblyname.txt
export families=${seqdb}/protein_families
export nrprotali=${protali}/nr_protfam_clustalo_alignments
export cdsalifastacodedir=${protali}/full_cdsfam_alignments_species_code
export protalifastacodedir=${protali}/full_protfam_alignments_species_code
export coretree=${coregenome}/raxml_tree
export mlgenetrees=${genetrees}/raxml_trees
export goterms=${funcannot}/GeneOntology
export claderefgodir=${goterms}/clade_go_term_reference_sets
export dirgoenrichcladespecore=${goterms}/clade_go_term_enriched_cladespecific_vs_coregenome
export dirgoenrichcladespepan=${goterms}/clade_go_term_enriched_cladespecific_vs_pangenome
#export dirgotablescladespe=${orthomatrad}_specific_genes.tables_byclade_goterms_pathways
# sub folders that depend on the gene tree clade collapsing option
export colalinexuscodedir=${genetrees}/${chaintype}_cdsfam_alignments_species_code
export recs=${alerec}/${chaintype}_${recmethod}_recs
export compoutdir=${comparerecs}/${parsedreccol}

# other secondary variables
export ngenomes=$(ls -A "${indata}/assemblies/" 2> /dev/null | wc -l)
if [ -d "${customassemb}" ] ; then
  export straininfo=${customassemb}/strain_infos_${ptgdbname}.txt
else
  if [ ! -z "${customstraininfo}" ] ; then
    export straininfo=${customstraininfo}
  else
    export straininfo=''
  fi
fi
export sqldbname=${ptgdbname,,}
export sqldb=${database}/${sqldbname}
export allfaarad=${seqdb}/all_proteomes
export mmseqsclout=${families}/$(basename ${allfaarad}.nr).mmseqs_clusterdb_default
export protfamseqs=${mmseqsclout}_clusters_fasta
export protorfanclust="${famprefix}P000000"
export cdsorfanclust="${famprefix}C000000"
if [ ! -z "$userreftree" ] ; then
  export coretreerad=${coregenome}/user-defined_reference_tree_${ptgdbname}
else
  export coretreerad=${coregenome}/core-genome-based_reference_tree_${ptgdbname}
fi
export nrbesttopo=${coretreerad}.topology
export nrbesttree=${coretreerad}.branlen
export nrbiparts=${coretreerad}.supports
export nrrootedtree=${coretreerad}.rooted
export speciestree=${coretreerad}.full
if [[ -z "${pseudocoremingenomes}" || "${pseudocoremingenomes}" == "${ngenomes}" ]] ; then
  export pseudocore='strict-core-unicopy'
else
  export pseudocore=pseudo-core-${pseudocoremingenomes}-unicopy
fi
export treename=${pseudocore}_concat_${coreseqtype}_${ngenomes}-genomes_${ptgdbname}
export pseudocorealn=${coregenome}/${treename}.aln
if [[ "${chaintype}" == 'fullgenetree' ]] ; then
  export collapsecond='nocollapse'
  export replmethod='noreplace'
#  export colmlgenetrees=${mlgenetrees}/rootedTree
  export colmlgenetrees=${mlgenetrees}/bipartitions 
  export coltreechains=${genetrees}/full_ML_genetrees
else
  export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}
  export replmethod='replaceCCinGasinS-collapsePOPinSnotinG'
  export colmlgenetrees=${colalinexuscodedir}/${collapsecond}/collapsed_ML_genetrees
  export coltreechains=${genetrees}/replaced_ML_genetrees
fi
export mboutputdir=${bayesgenetrees}/${collapsecond}
if [ -z "${pathtoipscan}" ] ; then
  export ipscanexe=interproscan
else
  export ipscanexe=${ptgdb}/interproscan
  rm -f ${ipscanexe}
  if [ -d "${pathtoipscan}" ] ; then
    ln -s ${pathtoipscan}/interproscan ${ipscanexe}
  else
	ln -s ${pathtoipscan} ${ipscanexe}
  fi
fi
export IPversion=$(${ipscanexe} --version 2> /dev/null | head -n 1 | sed -e 's/InterProScan version //')
if [ ! -z "${IPversion}" ] ; then
  export interpro=${funcannot}/InterProScan_${IPversion}
else
  export interpro=''
fi
