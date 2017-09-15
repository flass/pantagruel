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
export protali=$rapdb/02.clustalo_alignments
export database=$rapdb/03.database
export coregenome=$rapdb/04.core_genome

# sub folders
export complete=${assmetadata}/complete_genomes
export allcomplete=${complete}/complete_genomes
export faacomplete=$seqdb/all_complete_proteomes.faa
export nrfaacomplete=$seqdb/all_complete_proteomes.nr.faa
export families=$seqdb/protein_families
export nrprotali=$protali/nr_protfam_clustalo_alignments
export alifastacodedir=${protali}/full_cdsfam_alignments_species_code

# other variables
export sqldbname=${entdbname,,}
export protorfanclust="${famprefix}P000000"

# other variables conditonal on prior creation of files
if [ -e $nrfaacomplete ] ; then
 export mmseqsclout=${families}/$(basename $nrfaacomplete).mmseqs-clusterdb_default
 export protfamseqs=${mmseqsclout}_clusters_fasta
fi
if [ -e ${complete}/complete_genomes_metadata.tab ] ; then
 export ngenomes=$((`wc -l ${complete}/complete_genomes_metadata.tab | cut -d' ' -f1` - 1))
fi
