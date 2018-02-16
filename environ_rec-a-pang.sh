#!/bin/bash

## variables automaticly derived from user-defined environment variables $raproot and $rapdb
# head folders
export ncbiass=$raproot/NCBI/Assembly
export ncbitax=$raproot/NCBI/Taxonomy
export rapdb=$raproot/$rapdbname
export raplogs=$rapdb/logs
export raptmp=$rapdb/tmp
export indata=$rapdb/00.input_data
export seqdb=$rapdb/01.seqdb
export protali=$rapdb/02.gene_alignments
export database=$rapdb/03.database
export coregenome=$rapdb/04.core_genome
export genetrees=$rapdb/05.gene_trees


# sub folders
export complete=${assmetadata}/complete_genomes
export allcomplete=${complete}/complete_genomes
export faacomplete=$seqdb/all_complete_proteomes.faa
export nrfaacomplete=$seqdb/all_complete_proteomes.nr.faa
export families=$seqdb/protein_families
export nrprotali=$protali/nr_protfam_clustalo_alignments
export alifastacodedir=${protali}/full_cdsfam_alignments_species_code
export mlgenetrees=${genetrees}/raxml_trees
export colalinexuscodedir=${protali}/collapsed_cdsfam_alignments_species_code
export collapsed_genetrees=${genetrees}/collapsed_mrbayes_trees
export coltreechains=$alerec/collapsed_tree_chains
export colrecs=${alerec}/collapsed_recs

# other variables
export sqldbname=${entdbname,,}
export protorfanclust="${famprefix}P000000"
export cdsorfanclust="${famprefix}C000000"
export pseudocore=pseudo-core-${pseudocoremingenomes}-unicopy
export pseudocorealn=${coregenome}/${pseudocore}_concat_cds.aln
export coretree=${coregenome}/raxml_tree
export mainresulttag=rootedTree
cladesupp=70
subcladesupp=35
criterion='bs'
withinfun='median'
export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}
export colmethod='replaceCCinGasinS-collapsePOPinSnotinG'

# other variables conditonal on prior creation of files
if [ -e $nrfaacomplete ] ; then
 export mmseqsclout=${families}/$(basename $nrfaacomplete).mmseqs-clusterdb_default
 export protfamseqs=${mmseqsclout}_clusters_fasta
fi
if [ -e ${complete}/complete_genomes_metadata.tab ] ; then
 export ngenomes=$((`wc -l ${complete}/complete_genomes_metadata.tab | cut -d' ' -f1` - 1))
 export treename=${pseudocore}_concat_cds_${ngenomes}entero
 export nrspeciestree=$coretree/RAxML_bestTree.${treename}.MADrooted
 export nrspeciestreeBS=$coretree/RAxML_bipartitions.${treename}
fi


