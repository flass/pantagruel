#!/bin/bash

# logging variables and functions
alias dateprompt="date +'[%Y-%m-%d %H:%M:%S]'"
datepad="                     "

## variables to be automaticly replaced by user-defined environment variables
export raproot='REPLACEraproot'
export rapdbname='REPLACErapdbname'
export ptgscripts='REPLACEptgscripts'
export famprefix='REPLACEfamprefix'
export pseudocoremingenomes='REPLACEpseudocoremingenomes'

# head folders
export ncbiass=${raproot}/NCBI/Assembly
export ncbitax=${raproot}/NCBI/Taxonomy
export rapdb=${raproot}/${rapdbname}
export customassemb=${raproot}/user_genomes
export raplogs=${rapdb}/logs
export raptmp=${rapdb}/tmp
export indata=${rapdb}/00.input_data
export seqdb=${rapdb}/01.seqdb
export protali=${rapdb}/02.gene_alignments
export database=${rapdb}/03.database
export coregenome=${rapdb}/04.core_genome
export genetrees=${rapdb}/05.gene_trees


# sub folders
export annot=${customassemb}/prokka_annotation
export genomeinfo=${indata}/genome_infos
export assemblies=${indata}/assemblies
export protfamseqs=${mmseqsclout}_clusters_fasta
export families=${seqdb}/protein_families
export nrprotali=${protali}/nr_protfam_clustalo_alignments
export cdsalifastacodedir=${protali}/full_cdsfam_alignments_species_code
export protalifastacodedir=${protali}/full_protfam_alignments_species_code
export mlgenetrees=${genetrees}/raxml_trees
export colalinexuscodedir=${protali}/collapsed_cdsfam_alignments_species_code
export collapsed_genetrees=${genetrees}/collapsed_mrbayes_trees
export coltreechains=${alerec}/collapsed_tree_chains
export colrecs=${alerec}/collapsed_recs

# other variables
export straininfo=${customassemb}/strain_infos.txt
export sqldbname=${rapdbname,,}
export allfaarad=${seqdb}/all_proteomes
export mmseqsclout=${families}/$(basename ${allfaarad}.nr).mmseqs_clusterdb_default
export protorfanclust="${famprefix}P000000"
export cdsorfanclust="${famprefix}C000000"
export pseudocore=pseudo-core-${pseudocoremingenomes}-unicopy
export coretree=${coregenome}/raxml_tree
export mainresulttag=rootedTree
cladesupp=70
subcladesupp=35
criterion='bs'
withinfun='median'
export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}
export colmethod='replaceCCinGasinS-collapsePOPinSnotinG'

# other variables conditonal on prior creation of files
if [ -e ${complete}/complete_genomes_metadata.tab ] ; then
 export ngenomes=$((`wc -l ${genomeinfo}/metadata_${rapdbname}/metadata.tab | cut -d' ' -f1` - 1))
 export treename=${rapdbname}_${pseudocore}-concat-prot_${ngenomes}-genomes
 export pseudocorealn=${coregenome}/${treename}.aln
 export nrspeciestree=${coretree}/RAxML_bestTree.${treename}.MADrooted
 export nrspeciestreeBS=${coretree}/RAxML_bipartitions.${treename}
fi


