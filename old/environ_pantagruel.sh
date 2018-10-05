#!/bin/bash

# logging variables and functions
alias dateprompt="date +'[%Y-%m-%d %H:%M:%S]'"
datepad="                     "

## primary variables

# variables to be automaticly replaced by user-defined environment variables at init stage
export ptgroot='REPLACEptgroot'
export ptgdbname='REPLACEptgdbname'
export ptgscripts='REPLACEptgscripts'
export pseudocoremingenomes='REPLACEpseudocoremingenomes'
export famprefix='REPLACEfamprefix'


# parameters to be set; default values:
export chaintype='collapsed'
export cladesupp=70
export subcladesupp=35
export criterion='bs'
export withinfun='median'
export ALEalgo='ALEml_undated'
export recsamplesize=1000
export evtypeparse='ST'
export minevfreqparse=0.1
export minevfreqmatch=0.5
export minjoinevfreqmatch=1.0
export maxreftreeheight=0.25
export colmethod='replaceCCinGasinS-collapsePOPinSnotinG'
export mainresulttag=rootedTree

## secondary variables
# lib/module path
export PYTHONPATH=$PYTHONPATH:${ptgscripts}
# head folders
export ncbiass=${ptgroot}/NCBI/Assembly
export ncbitax=${ptgroot}/NCBI/Taxonomy
export ptgdb=${ptgroot}/${ptgdbname}
export customassemb=${ptgroot}/user_genomes
export ptglogs=${ptgdb}/logs
export ptgtmp=${ptgdb}/tmp
export indata=${ptgdb}/00.input_data
export seqdb=${ptgdb}/01.seqdb
export protali=${ptgdb}/02.gene_alignments
export database=${ptgdb}/03.database
export coregenome=${ptgdb}/04.core_genome
export genetrees=${ptgdb}/05.gene_trees
export alerec=${ptgdb}/06.ALE_reconciliation
export comparerecs=${entdb}/07.compare_scenarios
export orthogenes=${ptgdb}/08.orthologs

# sub folders
export annot=${customassemb}/prokka_annotation
export genomeinfo=${indata}/genome_infos
export assemblies=${indata}/assemblies
export protfamseqs=${mmseqsclout}_clusters_fasta
export families=${seqdb}/protein_families
export nrprotali=$protali/nr_protfam_clustalo_alignments
export cdsalifastacodedir=${protali}/full_cdsfam_alignments_species_code
export protalifastacodedir=${protali}/full_protfam_alignments_species_code
export mlgenetrees=${genetrees}/raxml_trees
# sub folders that depend on the gene tree clade collapsing option
export colalinexuscodedir=${protali}/${chaintype}_cdsfam_alignments_species_code
export bayesgenetrees=${genetrees}/${chaintype}_mrbayes_trees
export coltreechains=${alerec}/${chaintype}_tree_chains
export recs=${alerec}/${chaintype}_recs

# other secondary variables
export straininfo=${customassemb}/strain_infos.txt
export sqldbname=${ptgdbname,,}
export sqldb=${database}/${sqldbname}
export allfaarad=${seqdb}/all_proteomes
export mmseqsclout=${families}/$(basename ${allfaarad}.nr).mmseqs_clusterdb_default
export protorfanclust="${famprefix}P000000"
export cdsorfanclust="${famprefix}C000000"
export pseudocore=pseudo-core-${pseudocoremingenomes}-unicopy
export pseudocorealn=${coregenome}/${pseudocore}_concat_cds.aln
export coretree=${coregenome}/raxml_tree
export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}

# other secondary variables that depend on prior creation of files
if [ -e ${genomeinfo}/metadata_${rapdbname}/metadata.tab ] ; then
 export ngenomes=$((`wc -l ${genomeinfo}/metadata_${rapdbname}/metadata.tab | cut -d' ' -f1` - 1))
 export treename=${pseudocore}_concat_prot_${ngenomes}-genomes_${rapdbname}
 export pseudocorealn=${coregenome}/${treename}.aln
 nrbesttree=${coretree}/RAxML_bestTree.${treename}
 rootingmethod='outgroup'
 nrrootedtree=${nrbesttree}.${rootingmethod}rooted
 nrbiparts=${nrbesttree/bestTree/bipartitions}
 export speciestree=${nrrootedtree}.full
fi


