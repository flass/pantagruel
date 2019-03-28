#!/bin/bash

# logging variables and functions
promptdate () {
  echo $(date +'[%Y-%m-%d %H:%M:%S]') $1
}
export -f promptdate
datepad="                      "

## primary variables

# OPTIONALLY, values of primary variables below can be manually modified here;
# once edited, MAKE SURE TO SAVE THE FILE IN ANOTHER LOCATION
# and use it as init_file argument for command `pantagruel init init_file`
# NOTE this will prevent storing the value passed through options of `pantagruel [options] init` call.
# stored values can however be overridden by specifying the option when calling a specific task
# (only valid for task-relevant options, e.g. the use of -H option when calling 'genetrees' task will update the $hpcremoteptgrootvariable)

# variables to be automaticly set during the pantagruel run
# the first lot are derived from arguments passed as options in call `pantagruel [options] init`,
# for instance, options `-d TheDatabaseName` and `-r /where/it/is/built`
# will lead to setting variables ptgdbname='TheDatabaseName' and ptgroot='/where/it/is/built', respectively.
export ptgroot='REPLACEptgroot'                            # root folder where to build the database
export ptgdbname='REPLACEptgdbname'                        # name of dataabse
export ptgrepo='REPLACEptgrepo'                            # location (folder) of Pantagruel software
export myemail='REPLACEmyemail'                            # user identity (better use e-amil address)
export famprefix='REPLACEfamprefix'                        # gene family prefix
export ncbitax='REPLACEncbitax'                            # folder of up-to-date NCBI Taxonomy database
export ncbiass='REPLACEncbiass'                            # folder of RefSeq genomes to include in the study
export listncbiass='REPLACElistncbiass'                    # list of accessions of RefSeq genomes to include in the study
export customassemb='REPLACEcustomassemb'                  # folder of custom genome assemblies to include in the study
export refass='REPLACErefass'                              # folder of reference (RefSeq) genomes only to use as reference for the annotation of custom genome assemblies
export listrefass='REPLACElistrefass'                      # list of accessions of reference (RefSeq) genomes only to use as reference for the annotation of custom genome assemblies
export ngenomes=REPLACEngenomes                            # the count of genomes in the dataset
export coreseqtype='REPLACEcoreseqtype'                    # either 'cds' or 'protein'
export pseudocoremingenomes=REPLACEpseudocoremingenomes    # the minimum number of genomes in which a gene family should be present to be included in the pseudo-core genome gene set
export userreftree='REPLACEuserreftree'                    # possible user-provided reference tree
export poplgthresh='REPLACEpoplgthresh'                    # parameter to define populations of genomes in the reference tree (stem branch length threshold, default value depends on coreseqtype)
export poplgleafmul='REPLACEpoplgleafmul'                  # parameter to define populations of genomes in the reference tree (multiplier to the former in case it is a leaf, default 1.5)
export popbsthresh='REPLACEpopbsthresh'                    # parameter to define populations of genomes in the reference tree (stem branch support threshold, default 80)
export chaintype='REPLACEchaintype'                        # whether gene trees will be collapsed ('collapsed', if -c option enabled) or not ('fullgenetree', default)
# non-default parameters for gene trees collapsing derived from -C option value (passed to init script via ${collapseCladeParams}): 
export cladesupp=REPLACEcladesupp                          # - clade criterion trheshold (int)
export subcladesupp=REPLACEsubcladesupp                    # - wihtin-clade criterion trheshold (int)
export criterion='REPLACEcriterion'                        # - criterion (branch support: 'bs', branch length 'lg')
export withinfun='REPLACEwithinfun'                        # - aggregate function for testing within the clade ('min', 'max', 'mean', 'median')
export hpcremoteptgroot='REPLACEhpcremoteptgroot'          # if not empty nor 'none', will use this server address to send data and scripts to run heavy computions there 
# the rest have fixed values that can be modified here or overriden with certain task-specific options
# default values are:
export cladesuppdef=70
export subcladesuppdef=35
export criteriondef='bs'
export withinfundef='median'
export ALEalgo='ALEml_undated'
export recsamplesize=1000
export evtypeparse='ST'
export minevfreqparse=0.1
export minevfreqmatch=0.5
export minjoinevfreqmatch=1.0
export maxreftreeheight=0.25
export colmethod='replaceCCinGasinS-collapsePOPinSnotinG'
export ncorebootstrap=200
export mainresulttag='rootedTree'
export rootingmethod='treebalance'
# only relevant if user-defined genomes are provided
export assembler="somesoftware"
export seqcentre="somewhere"
export refgenus="Reference"

### BEWARE!!! value of variables below is unsafe to modify

## secondary variables
# lib/module path
export ptgscripts=${ptgrepo}/scripts
export PYTHONPATH=$PYTHONPATH:${ptgrepo}/python_libs
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
export alerec=${ptgdb}/07.ALE_reconciliation
export orthogenes=${ptgdb}/08.orthologs
export comparerecs=${entdb}/09.compare_scenarios

# sub folders
export contigs=${customassemb}/contigs
export annot=${customassemb}/annotation
export gblikeass=${customassemb}/genbank-format_assemblies
export genomeinfo=${indata}/genome_infos
export assemblies=${indata}/assemblies
export families=${seqdb}/protein_families
export nrprotali=$protali/nr_protfam_clustalo_alignments
export cdsalifastacodedir=${protali}/full_cdsfam_alignments_species_code
export protalifastacodedir=${protali}/full_protfam_alignments_species_code
export coretree=${coregenome}/raxml_tree
export mlgenetrees=${genetrees}/raxml_trees
# sub folders that depend on the gene tree clade collapsing option
export colalinexuscodedir=${genetrees}/${chaintype}_cdsfam_alignments_species_code
export bayesgenetrees=${genetrees}/${chaintype}_mrbayes_trees
export coltreechains=${genetrees}/${chaintype}_tree_chains
export recs=${alerec}/${chaintype}_recs
export goterms=${funcannot}/GeneOntology
export claderefgodir=${goterms}/clade_go_term_reference_sets
export dirgotablescladespe=${orthomatrad}_specific_genes.tables_byclade_goterms_pathways
export dirgoenrichcladespecore=${goterms}/clade_go_term_enriched_cladespecific_vs_coregenome
export dirgoenrichcladespepan=${goterms}/clade_go_term_enriched_cladespecific_vs_pangenome
export compoutdir=${comparerecs}/${parsedreccol}

# other secondary variables
export ngenomes=$(ls -A "${indata}/assemblies/" 2> /dev/null | wc -l)
export straininfo=${customassemb}/strain_infos_${ptgdbname}.txt
export sqldbname=${ptgdbname,,}
export sqldb=${database}/${sqldbname}
export allfaarad=${seqdb}/all_proteomes
export mmseqsclout=${families}/$(basename ${allfaarad}.nr).mmseqs_clusterdb_default
export protfamseqs=${mmseqsclout}_clusters_fasta
export protorfanclust="${famprefix}P000000"
export cdsorfanclust="${famprefix}C000000"
if [ "$userreftree" != "REPLACEuserreftree" ] ; then
  export coretreerad=${coregenome}/user-defined_reference_tree_${ptgdbname}
else
  export coretreerad=${coregenome}/core-genome-based_reference_tree_${ptgdbname}
fi
export nrbesttree=${coretreerad}.topology
export nrbiparts=${coretreerad}.supports
export nrrootedtree=${coretreerad}.rooted
export speciestree=${coretreerad}.full
if [[ "${pseudocoremingenomes}" == "${ngenomes}" ]] ; then
  export pseudocore='strict-core-unicopy'
else
  export pseudocore=pseudo-core-${pseudocoremingenomes}-unicopy
fi
export treename=${pseudocore}_concat_${coreseqtype}_${ngenomes}-genomes_${ptgdbname}
export pseudocorealn=${coregenome}/${treename}.aln
export collapsecond=${criterion}_stem${cladesupp}_within${withinfun}${subcladesupp}
export IPversion=$(interproscan --version | head -n 1 | sed -e 's/InterProScan version //')
export interpro=${funcannot}/InterProScan_${IPversion}


