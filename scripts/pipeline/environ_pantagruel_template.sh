#!/bin/bash

# logging variables and functions
alias dateprompt="date +'[%Y-%m-%d %H:%M:%S]'"
datepad="                     "

## primary variables

# variables to be automaticly replaced by user-defined environment variables
# at init stage
export ptgroot='REPLACEptgroot'
export ptgdbname='REPLACEptgdbname'
export ptgrepo='REPLACEptgrepo'
export myemail='REPLACEmyemail'
export famprefix='REPLACEfamprefix'
# at 00.input_data stage
export ngenomes='REPLACEngenomes'
# at 04.core_genome stage
export pseudocoremingenomes='REPLACEpseudocoremingenomes'

# OPTIONALLY, values of parameters below can be manually modified here;
# once edited, MAKE SURE TO SAVE THE FILE IN ANOTHER LOCATION
# and use it as init_file argument for command `pantagruel init init_file`
# default values are:
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
export mainresulttag='rootedTree'
export rootingmethod='treebalance'
# only relevant if user-defined genomes are provided
export assembler="somemachine"
export seqcentre="somewhere"

### BEWARE!!! definition of variables below is unsafe to modify

## secondary variables
# lib/module path
export ptgscripts=${ptgrepo}/scripts
export PYTHONPATH=$PYTHONPATH:${ptgrepo}/python_libs
# head folders
export ptgdb=${ptgroot}/${ptgdbname}
export ncbiass=${ptgroot}/NCBI/Assembly
export ncbitax=${ptgroot}/NCBI/Taxonomy
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
export families=${seqdb}/protein_families
export nrprotali=$protali/nr_protfam_clustalo_alignments
export cdsalifastacodedir=${protali}/full_cdsfam_alignments_species_code
export protalifastacodedir=${protali}/full_protfam_alignments_species_code
export coretree=${coregenome}/raxml_tree
export mlgenetrees=${genetrees}/raxml_trees
# sub folders that depend on the gene tree clade collapsing option
export colalinexuscodedir=${protali}/${chaintype}_cdsfam_alignments_species_code
export bayesgenetrees=${genetrees}/${chaintype}_mrbayes_trees
export coltreechains=${alerec}/${chaintype}_tree_chains
export recs=${alerec}/${chaintype}_recs

# other secondary variables
export straininfo=${customassemb}/${ptgdbname}_strain_infos.txt
export sqldbname=${ptgdbname,,}
export sqldb=${database}/${sqldbname}
export allfaarad=${seqdb}/all_proteomes
export mmseqsclout=${families}/$(basename ${allfaarad}.nr).mmseqs_clusterdb_default
export protfamseqs=${mmseqsclout}_clusters_fasta
export protorfanclust="${famprefix}P000000"
export cdsorfanclust="${famprefix}C000000"
export pseudocore=pseudo-core-${pseudocoremingenomes}-unicopy
export treename=${pseudocore}_concat_prot_${ngenomes}-genomes_${ptgdbname}
export pseudocorealn=${coregenome}/${treename}.aln
export nrbesttree=${coretree}/RAxML_bestTree.${treename}
export nrbiparts=${nrbesttree/bestTree/bipartitions}
export nrrootedtree=${nrbiparts}.${rootingmethod}rooted
export speciestree=${nrrootedtree}.full
export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}



