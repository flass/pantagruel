#!/bin/bash
# head folders
export ncbiass=$raproot/NCBI/Assembly
export ncbitax=$raproot/NCBI/Taxonomy
export rapdb=$raproot/$rapdbname
export raplogs=$rapdb/logs
export raptmp=$rapdb/tmp
export indata=$rapdb/00.input_data
export seqdb=$rapdb/01.seqdb
export protali=$rapdb/02.clustalo_alignments
# sub folders
export complete=${assmetadata}/complete_genomes
export allcomplete=${complete}/complete_genomes
export faacomplete=$seqdb/all_complete_proteomes.faa
export nrfaacomplete=$seqdb/all_complete_proteomes.nr.faa
export families=$seqdb/protein_families
mmseqsclout=${families}/$(basename $nrfaacomplete).mmseqs-clusterdb_default
export protfamseqs=${mmseqsclout}_clusters_fasta
nrprotali=$protali/nr_protfam_clustalo_alignments
