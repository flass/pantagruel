#!/bin/bash

# This is the template for generating a Pantagruel database configuration file

# Variables with values starting with REPLACE are to be automaticly set during the pantagruel intit run
# based on arguments passed as options in call `pantagruel [options] init`,
# for instance, options `-d TheDatabaseName` and `-r /where/it/is/built`
# will lead to setting variables ptgdbname='TheDatabaseName' and ptgroot='/where/it/is/built', respectively.
#
# OPTIONALLY, values of primary variables below can be manually modified here;
# once edited, MAKE SURE TO SAVE THE FILE IN ANOTHER LOCATION THAN THE PANTAGRUEL SOFTWARE REPOSITORY TEMPLATE FILE
# and use it as init_file argument for command `pantagruel init init_file`
# NOTE this will prevent storing the value passed through options of `pantagruel [options] init` call.
# stored values can however be overridden by specifying the option when calling a specific task
# (only valid for task-relevant options, e.g. the use of -H option when calling 'genetrees' task will update the $hpcremoteptgrootvariable)

# the init command with which the config file was created
ptginitcmd='REPLACEptginitcmd'

# location (folder) of Pantagruel software that was used
export ptgrepo='REPLACEptgrepo'                            
# derive paths to Pantagruel scripts and Python modules
export ptgscripts=${ptgrepo}/scripts
export PYTHONPATH=${ptgrepo}/python_libs
# database parameters (primary variables)
export ptgroot='REPLACEptgroot'                            # root folder where to build the database
export ptgdbname='REPLACEptgdbname'                        # name of dataabse
export ptgversinit='REPLACEptgversinit'                    # current version of Pantagruel software
export myemail='REPLACEmyemail'                            # user identity (better use e-amil address)
export famprefix='REPLACEfamprefix'                        # gene family prefix
export ncbitax='REPLACEncbitax'                            # folder of up-to-date NCBI Taxonomy database
export ncbiass='REPLACEncbiass'                            # folder of RefSeq genomes to include in the study
export listncbiass='REPLACElistncbiass'                    # list of accessions of RefSeq genomes to include in the study
export customassemb='REPLACEcustomassemb'                  # folder of custom genome assemblies to include in the study
export refass='REPLACErefass'                              # folder of reference (RefSeq) genomes only to use as reference for the annotation of custom genome assemblies
export listrefass='REPLACElistrefass'                      # list of accessions of reference (RefSeq) genomes only to use as reference for the annotation of custom genome assemblies
export coreseqtype='REPLACEcoreseqtype'                    # either 'cds' or 'protein'
export pseudocoremingenomes=REPLACEpseudocoremingenomes    # the minimum number of genomes in which a gene family should be present to be included in the pseudo-core genome gene set
export userreftree='REPLACEuserreftree'                    # possible user-provided reference tree
export poplgthresh='REPLACEpoplgthresh'                    # parameter to define populations of genomes in the reference tree (stem branch length threshold, default value depends on coreseqtype)
export poplgleafmul='REPLACEpoplgleafmul'                  # parameter to define populations of genomes in the reference tree (multiplier to the former in case it is a leaf, default 1.5)
export popbsthresh='REPLACEpopbsthresh'                    # parameter to define populations of genomes in the reference tree (stem branch support threshold, default 80)
export rootingmethod='REPLACErootingmethod'                # rooting method for core-genome tree
export snpali='REPLACEsnpali'                              # restrict core-genome alignment to SNPs
export chaintype='REPLACEchaintype'                        # whether gene trees will be collapsed ('collapsed', if -c option enabled) or not ('fullgenetree', default)
export genefamlist='REPLACEgenefamlist'                    # list of gene families for which computation of gene trees and all subsequent analyses will be restricted (default: no restriction)
export recmethod='REPLACErecmethod'                        # genetree/species tree reconciliation method: 'ALE' or 'ecceTERA'
# non-default parameters for gene trees collapsing derived from -C option value (passed to init script via ${collapseCladeParams}): 
export cladesupp=REPLACEcladesupp                          # - clade criterion trheshold (int)
export subcladesupp=REPLACEsubcladesupp                    # - wihtin-clade criterion trheshold (int)
export criterion='REPLACEcriterion'                        # - criterion (branch support: 'bs', branch length 'lg')
export withinfun='REPLACEwithinfun'                        # - aggregate function for testing within the clade ('min', 'max', 'mean', 'median')
export hpcremoteptgroot='REPLACEhpcremoteptgroot'          # if not empty nor 'none', will use this server address to send data and scripts to run heavy computions there 
export maxreftreeheight='REPLACEmaxreftreeheight'          # restict events younger than that age (comprised in [0.0; 1.0]) on the species tree for gene co-evolution scoring
export updatedbfrom='REPLACEupdatedbfrom'                  # the current pantagruel database is an update from that found at this path
export customstraininfo='REPLACEcustomstraininfo'          # optional custom strain information file
export pathtoipscan='REPLACEpathtoipscan'                  # optional path to interproscan executable
## other parameters have default values defined in the generic source file environ_pantagruel_defaults.sh
source ${ptgscripts}/pipeline/environ_pantagruel_defaults.sh
## these defalts can be overriden by uncommenting the relevant line below and editing the variable's value
## or (recomended for changes to last past calls to `pantagruel --refresh init`):
## create a file '${ptgroot}/${ptgdbname}/user_environ_pantagruel_${ptgdbname}.sh' containing the `export variable=value` commands
# default values are:
# Prokka annotation parameters (only relevant if custom genome assemblies are provided):
#~ export assembler="somesoftware"
#~ export seqcentre="somewhere"
#~ export refgenus="Reference"
# species tree inference parameters
#~ export ncorebootstrap=200
# gene tree inference parameters
#~ export mainresulttag='rootedTree'
# gene trees collapsing DEFAULT values (used when -C option is NOT present in init call)
#~ export cladesuppdef=70
#~ export subcladesuppdef=35
#~ export criteriondef='bs'
#~ export withinfundef='median'
# gene tree/species tree reconciliation inference parameters
#~ export ALEalgo='ALEml'
#~ export ecceTERAalgo='amalgamate'
#~ export recsamplesize=1000
# gene tree/species tree reconciliation parsing parameters for co-evolution analysis
#~ export evtypeparse='ST'
#~ export minevfreqparse=0.1
#~ export minevfreqmatch=0.5
#~ export minjoinevfreqmatch=1.0
#~ export maxreftreeheight=0.25
userparams="${ptgroot}/${ptgdbname}/user_environ_pantagruel_${ptgdbname}.sh"
if [ -s "${userparams}" ] ; then
  echo "Warning: will use user-defined values for Pantagruel environment variables, as deined in '${userparams}':"
  cat ${userparams}
  source ${userparams}
fi

# secondary vars are defined based on the above
source ${ptgscripts}/pipeline/environ_pantagruel_secondaryvars.sh
# load shared functions
source ${ptgscripts}/pipeline/pantagruel_pipeline_functions.sh
