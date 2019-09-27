#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################
# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 01 October 2018.

ptgcmdargs=$(echo ${@} | sed -e "s/'/\'/g")
locexec=$(readlink -e $0)
export ptgrepo=${locexec%%/scripts/pipeline/*}
cd ${ptgrepo} ; export ptgversion=$(git log | head -n 1 | awk '{ print $2 }') ; cd - > /dev/null
echo "This is Pantagruel pipeline version ${ptgversion} using source code from repository '$ptgrepo'"
export ptgscripts=${ptgrepo}/scripts

usage (){
  echo "Usage: pantagruel -d db_name -r root_dir {-A refseq_assembly_folder | -a custom_assembly_folder} [other init options] init"
  echo " or"
  echo "Usage: pantagruel -i config_file TASK1 [TASK2, [...]]"
}

usagelong (){
  usage
  echo "# for Pantagruel task 0-9:"
  echo ""
  echo "_only one mandatory option_"
  echo ""
  echo "    -i|--initfile     Pantagruel configuration file"
  echo "                        this file is generated at init stage, from the specified options."
  echo ""
  echo "  _facultative runtime options_"
  echo ""
  echo "    -F|--FORCE        FORCE mode: will erase any pre-existing main folder for the task"
  echo "                        (default: off, pre-exisitance of a folder will result in an early error)"
  echo ""
  echo "    -R|--resume       try and resume the task from previous run that was interupted"
  echo "                        (for the moment only available for tasks: 04-07 i.e. 'functional', 'core', 'genetrees' and 'reconciliations')"
  echo ""
  echo "    -N|--threads      specify the number of threads to use for (some) parrallelizable tasks (defaults to the maximum available: $(nproc))"
  echo "                        (for the moment only available for tasks: 00,04,05 i.e. 'fetch', 'functional', 'core')"
  echo ""
  echo "    -z|--compress     will try and compress result file on the go (especially bulky files that won't be used much later"
  echo "                        This will induce possible extra decompressing/re-generating data computing time"
  echo "                        when resuming a task run with -R; avoid using compression when likely to have to resume later"
  echo "                        (for the moment only available for tasks: 00 i.e. 'fetch')"
  echo ""
  echo "# for Pantagruel task init:"
  echo ""
  echo "  _mandatory options_"
  echo ""
  echo "    -d|--dbname       database name"
  echo ""
  echo "    -r|--rootdir      root directory where to create the database; defaults to current folder"
  echo ""
  echo "    It is also necessary to specify an input genome dataset!"
  echo "    This is possible via -a, -A or -L options, or a mixture of them."
  echo ""
  echo "  _facultative options_"
  echo ""
  echo "    -i|--initfile     Pantagruel configuration file"
  echo "                        a file can be derived (i.e. manualy curated) from 'environment_pantagruel_template.sh' template."
  echo "                        Parameters values specified in this file will override other options"
  echo "                        Can also be combined alone with --refresh to update the software version used for an existing database."
  echo ""
  echo "    --refresh         (no value) Use in combination with the -i option above to simply refresh the configuration file"
  echo "                        (e.g. after an update of the software) the program will simply re-run the \`pantagruel [options] init\` command"
  echo "                        that has been previously used to generate the config file; hence there is no need to repeat any other option."
  echo "                        (even -d and - r options can be omitted if -i --refresh is used)"
  echo "                        Note that when options had quoted string arguments, unpredictable behaviour might occur;"
  echo "                        please verify the outcome in the regenerated config file."
  echo ""
  echo " General configuration options:"
  echo ""
  echo "    -I|--iam          database creator identity (e-mail address is preferred)"
  echo ""
  echo "    -f|--famprefix    alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; defaults to 'PANTAG'"
  echo "                        the chosen prefix will be appended with a 'P' for protein families and a 'C' for CDS families."
  echo ""
  echo " Input options:"
  echo ""
  echo "    -T|--taxonomy     path to folder of taxonomy database flat files; defaults to \$rootdir/NCBI/Taxonomy_YYYY-MM-DD (suffix is today's date)"
  echo "                        if this is not containing the expected file, triggers downloading the daily dump from NCBI Taxonomy at task 00"
  echo ""
  echo "    -A|--refseq_ass   path to folder of source genome assembly __folders__, each containing flat files formated like NCBI Assembly RefSeq."
  echo "                        (no default value)"
  echo "                        The assembly folders (one per genome) can be obtained on https://www.ncbi.nlm.nih.gov/assembly,"
  echo "                        by making a keyword search to select a set of assemblies and downloading results with options:"
  echo "                          Source Database = 'RefSeq' and File type = 'All file types (including assembly-structure directory)'."
  echo "                        A simple archive 'genome_assemblies.tar' (as obtained from the NCBI website) can be placed in that folder."
  echo "                        If user genomes are also provided, these RefSeq assemblies will be used as reference for their annotation."
  echo ""
  echo "    --refseq_ass4annot folder path. Same principle as -A, but WILL NOT be used in the study, only as a reference to annotate user genomes"
  echo "                        (defaults to combined value of -A and -L options)"
  echo ""
  echo "    -L|--refseq_list  file path. Same as -A|--refseq_ass, but just provide a list of NCBI Assembly accession ids (file with one accession per row)"
  echo "                        Accession ids are of the form GCx_yyyyyyyyy.z with x = {A|F} for GenBank and RefSeq, respectively, and y an z are any digit."
  echo "                        These accessions will be fetched from the NCBI FTP site using lftp."
  echo "                        Note the LAST version of the accession will be returned, i.e. the trailing '.z' part of the accession id is ignored."
  echo "                        These assemblies saved into a folder named after the value of the option:"
  echo "                          for instance, \`-L /path/to/assemblist\` will save assembly folders in /path/to/assemblist_assemblies_from_ftp/."
  echo ""
  echo "    --refseq_list4annot file path. Same principle as -L, but WILL NOT be used in the study, only as a reference to annotate user genomes"
  echo "                        (defaults to the combined value of -A and -L options)"
  echo ""
  echo "    -a|--custom_ass   path to folder of user-provided genomes (no default value). The specified folder must contain:"
  echo "                       _mandatory_ "
  echo "                        - a 'contigs/' folder, where are stored multi-FASTA files of genome assemblies (one file per genome,"
  echo "                           with extension '.fa', '.fasta' or '.fas' ...). Fasta file names will be truncated by removing"
  echo "                           the '.fa' string and everything occuring after) and will be retained as the assembly_id (beware "
  echo "                           of names redundant with RefSeq assemblies)."
  echo "                        - a 'strain_infos_\$\{databasename\}.txt' TAB-delimited file describing the organism, with \$\{databasename\} the value of option -d"
  echo "                           columns should be headed with these fields (replace quotes and semicolons by tabs!):"
  echo "                            'assembly_id'; 'genus'; 'species'; 'strain'; 'taxid'; 'locus_tag_prefix'"
  echo "                          'assembly_id' must match the name of a contig file (e.g. 'seqProjID.fasta')"
  echo "                          'locus_tag_prefix' must match the prefix of ids given to CDS, proteins and genome regions (contigs)"
  echo "                          in potentially provided annotation files (see below)."
  echo "                          Note that for ensuring compatibility with dependencies (namely BioPython Nexus alignment parser"
  echo "                          and ALE reconciliation program), the characters '-' and '_' are forbiden in the 'locus_tag_prefix' field."
  echo "                       _optional_ "
  echo "                        - an 'annotation/' folder, where are stored annotation files: "
  echo "                          - one mandatory in GFF 3.0 file format (with a '.gff' extension);"
  echo "                           and optionally, the following files (with consistent ids!!):"
  echo "                          - one in GenBank flat file format (with a '.gbk' extension);"
  echo "                          - one in Fasta format containing CDS sequences (with a '.ffn' extension)."
  echo "                          - one in Fasta format containing matching protein sequences (with a '.faa' extension)."
  echo "                          These four files are produced when using Prokka for annotation; if at least one of the .gbk, .ffn or .faa"
  echo "                          are missing, all three will be derived from the .gff source. Each genome annotation file set must be stored"
  echo "                          in a separate folder, which name must match a contig file (e.g. 'seqProjID/' for 'seqProjID.fasta')."
  echo "                       NOTE: to ensure proper parsing, it is strongly advised that any provided annotation was generated with Prokka"
  echo "                       NOTE: to ensure uniform annotation of the dataset, it is advised to let Pantagruel annotate the contigs (calling Prokka)"
  echo ""
  echo "    -V|--env_var    (preferably double-quoted) string of the form: \"variable1=value1[,variable2=value2[,...]]'.\""
  echo "                     Will add these variables to the configuration file so they can be exported to the environment during tasks."
  echo "                     Can be useful to define custom values of generic variables, e.g. \"refgenus=Escherichia,seqcentre=Sanger_Institute\""
  echo ""
  echo " Output: core genome / reference phylogeny options:"
  echo ""
  echo "    -s|--pseudocore   integer, float <=1.0 or string. The minimum number or fraction of genomes in which a gene family should be present"
  echo "                        to be included in the pseudo-core genome, i.e. the gene set which alignments will be concatenated for reference tree search."
  echo "                        A non-numeric value will trigger an INTERACTIVE prompt for search of an optimal value at the begining of task 'core'."
  echo "                        Defaults to the total number of genomes (strict core genome set)."
  echo ""
  echo "    -t|--reftree      Newick format tree file path. Specifies a reference tree for reconciliation and clade-specific gene analyses;"
  echo "                        cancels the computation of tree from the concatenate of (pseudo-)core genome gene during task 'core'."
  echo ""
  echo "    --core_seqtype    {cds|prot} define the type of sequence that will be used to compute the (pseudo-)core genome tree (default to 'cds')"
  echo ""
  echo "    --pop_lg_thresh   real. Defines the threshold of branch length for delinating populations in the reference tree "
  echo "                        (default: 0.0005 for nucleotide alignment-based tree; 0.0002 for protein-based)"
  echo ""
  echo "    --pop_bs_thresh   real. Defines the threshold of branch support for delinating populations in the reference tree (default: 80)"
  echo ""
  echo "    --rooting         string. Defines the method to root the reference tree during task 5|core_genome_ref_tree. "
  echo "                        Possible values are 'treebalance', 'MAD' and 'outgroup:SPECIESCODELIST' (default: 'treebalance'),"
  echo "                        - 'treebalance' uses the '-f I' algorthm of RAxML to root the tree towards an optimal balance of branch lengths"
  echo "                           on either sides of the root;"
  echo "                        - 'MAD' uses the minimal ancestor deviation method described in \"Tria, et al. (2017) Nat. Ecol. Evol. 1, 0193\"."
  echo "                        - 'outgroup:SPECIESCODELIST' will root according tothe specified outgroup(s), with SPECIESCODELIST a comma-sperated list of species ids:"
  echo "                            'outgroup:SPECIESCODE' for rooting with a single species"
  echo "                            'outgroup:SPECIESCODE1,SPECIESCODE2,... for mutilple species (in which case their MRCA in  the tree will be the outgroup)"
  echo "                           Species ids can be either valid genome assembly ids of the relevant input genomes (typically a NCBI Assembly accession id),"
  echo "                           or internal genome identifiers that are specifically in the Pantagruel database but often match the relevant Uniprot organism code."
  echo "                           The mapping between genome accession ids and organism codes is given in the file '03.database/genome_codes.tab' generated during task 3."
  echo "                           To use codes, you may thus want run task 3 first, then run task init again with this option to regenerate the config file with "
  echo "                           the desired outgroup organism codes and only then run task 5."
  echo ""
  echo " Output: gene trees / reconciliations options:"
  echo ""
  echo "    -H|--submit_hpc   full address (hostname:/folder/location) of a folder on a remote high-performance computating (HPC) cluster server"
  echo "                        This indicate that computationally intensive tasks, including building the gene tree collection"
  echo "                        ('genetrees') and reconciling gene tree with species tree ('reconciliations') will be run"
  echo "                        on a HPC server (only Torque/PBS job submission system is supported so far)."
  echo "                        [support for core genome tree building ('core') remains to be implemented]."
  echo "                        Instead of running the computations, scripts for cluster job submission will be generated automatically."
  echo "                        Data and scripts will be transfered to the specified address (the database folder structure"
  echo "                        will be duplicated there, but only relevant files will be synced). Note that job submission"
  echo "                        scripts will need to be executed manually on the cluster server."
  echo "                        If set at init stage, this option will be maintained for all tasks. However, the remote address"
  echo "                        can be updated when calling a specific task; string 'none' cancels the HPC behaviour."
  echo ""
  echo "    -c|--collapse     (no value) enable collapsing the rake clades in the gene trees (strongly recomended in datasets of size > 50 genomes)."
  echo ""
  echo "    --collapse_param  quoted string. specify parameters for collapsing the rake clades in the gene trees."
  echo "                        A single-quoted, semicolon-delimited string containing variable definitions must be provided."
  echo "                        Default is equivalent to providing the following string:"
  echo "                           'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'"
  echo ""
  echo "    -g|--genefam_list Path to gene family list file. Resticts the computation of gene trees and all subsequent analyses to a list of gene families."
  echo "                        This impacts all task from 06 and forward. The list has to be one gene family identifier per line."
  echo "                        Gene family ids have to refer to existing ones in the database, and therefore can only be defined after the running of task 02."
  echo "                        It is therefore advised to first run the pipeline up to task 02 (or equally up to 05) without this option,"
  echo "                        and then to to set this paramter for the downstream computations."
  echo "                        This can be done by editing the value of 'genefamlist' variable in the configuration file or by using:"
  echo "                          pantagruel -i configfile --refresh -g genelist init  (note it is important that -g option be placed after the --refresh option)"
  echo "                        Reverting to the exhaustive computation behavior can be done similarly by setting 'genefamlist' variable to an empty value or by using:"
  echo "                          pantagruel -i configfile --refresh -g '' init"
  echo ""
  echo "# for any Pantagruel command calls:"
  echo ""
  echo "    -h|--help          print this help message and exit."
  echo ""
  echo "TASKs are to be picked among the following (equivalent digit/number/keywords are separated by a '|'):"
  echo "  init"
  echo "       database initiation: define shared environment variables."
  echo "         If followed by a file path, all env variables are set as in the specified source shell file"
  echo "  all"
  echo "       run all pipeline tasks (apart from database initiation), from fetching data to co-evolution inference."
  echo "       All tasks will be performed using parameters defined at the __init__ step, and default parameters otherwise."
  echo "  0|00|fetch|fetch_data"  
  echo "       fetch public genome data from NCBI sequence databases and annotate private genomes"
  echo "       also run a quick genome-to-genome distance estimation with MASH"
  echo "  1|01|homologous|homologous_seq_families"
  echo "       classify protein sequences into homologous families"
  echo "  2|02|align|align_homologous_seq"  
  echo "       align homologous protein sequences and translate alignemnts into coding sequences"
  echo "  3|03|sqldb|create_sqlite_db"    
  echo "       initiate SQL database and load genomic object relationships"  
  echo "  4|04|functional|functional_annotations"
  echo "       use InterProScan to functionally annotate proteins in the database, including with Gene Ontology and metabolic pathway terms"
  echo "  5|05|core|core_genome_ref_tree"
  echo "       select core-genome markers and compute reference tree"
  echo "  6|06|genetrees|gene_trees"
  echo "       compute gene tree"
  echo "  7|07|reconciliations"
  echo "       compute species tree/gene tree reconciliations"
  echo "  8|08|specific|clade_specific_genes"
  echo "       classify genes into orthologous groups (OGs) and search clade-specific OGs"
  echo "  9|09|coevolution"
  echo "       quantify gene co-evolution and build gene association network"
  echo "_________________________________________________________________________________________________________________"
  echo "     "
  echo "     Please note that Pantagruel is still under development and may evolve rapidly to fix bugs and solve issues. "
  echo "     It is strongly recomended to update this software regularly to integrate those fixes using \`git pull\` in the"
  echo "     pantagruel/ git repository folder."
  echo "     "
  echo "     If Pantagruel is updated in-between the running of tasks of a database project, it is higly recommended to "
  echo "     run the init task again (with the same options) before proceeding to the next steps, to ensure the environment "
  echo "     defined by the configuration file is compatible with the software. "
  echo "     (Note that regenerating the configuration file will not erase other data.)"
  echo "     A simple way to regenerate the configuration file under the same parameters is to use the command:"
  echo "        pantagruel -i previous_config_file --refresh init "
  echo "_________________________________________________________________________________________________________________"
  echo ""
}

ptgenvsetdefaults (){
    # Default values:
    export ptgscripts=${ptgrepo}/scripts
    if [ -z "$runmode" ] ; then
      export runmode="normal"
      echo "Default: set runmode to '$runmode'"
    fi
    if [ -z "$myemail" ] ; then
      export myemail="undisclosed"
      echo "Default: set identity to '$myemail'"
    fi
    if [ -z "$famprefix" ] ; then
      export famprefix="PANTAG"
      echo "Default: set gene family prefix to '$famprefix'"
    fi
    if [ -z "$ncbitax" ] ; then
     export ncbitax=${ptgroot}/NCBI/Taxonomy_$(date +'%Y-%m-%d')
     echo "Default: set NCBI Taxonomy source folder to '$ncbitax'"
    fi
    if [ -z "$chaintype" ] ; then
      export chaintype='fullgenetree'
      echo "Default: set gene tree type to '$chaintype'"
    fi
    if [ -z "$pseudocoremingenomes" ] ; then
     echo "Default: will use a strict core-genome gene set, i.e. genes present in a single copy in all the studied genomes."
     echo ""
     echo "!!! WARNING: strict core-genome definition can be very resctrictive, especially when including draft genome in the study."
     echo "You might prefer to use a pseudo-core genome definition instead, i.e. selecting gene present in a minimum fraction of genomes, for instance 98%."
     echo "A sensible threshold should avoid that selected genes have an approximately homogeneous distribution,"
     echo "notably that the absent fraction is not restricted to a few genomes. This threshold will thus depend on the dataset."
     echo "To choose a sensible value, AFTER TASK 03, you can run the INTERACTIVE script:"
     echo "'$ptgscripts/choose_min_genome_occurrence_pseudocore_genes.sh'"
     echo "and then manualy edit the value of variable 'pseudocoremingenomes' in the pantagruel configuration file."
    else
      if [[ "$pseudocoremingenomes" =~ ^[0-9]+$ ]]; then
      	 echo "'pseudocoremingenomes' variable is set to ${pseudocoremingenomes}; this integer value is interpreted as a number of genomes"
      elif [[ "$pseudocoremingenomes" =~ ^0*\.[0-9]+$ ]]; then
      	 echo "'pseudocoremingenomes' variable is set to ${pseudocoremingenomes}; this float value is interpreted as a fraction of total number of genomes"
      else		
      	 echo "'pseudocoremingenomes' variable is not set to a numeric value: '${pseudocoremingenomes}'; will run INTERACTIVE script to pick an appropriate one"
      fi
    fi
    if [ -z "$coreseqtype" ] ; then
     export coreseqtype='cds'
     echo "Default (only relevant to 'core' task): core sequence type is set to '$coreseqtype'"
    fi
    if [ -z "$poplgthresh" ] ; then
     export poplgthresh='default'
     echo "Default (only relevant to 'core' task): set population delination branch support threshold to $poplgthresh"
    fi
    if [ -z "$poplgleafmul" ] ; then
     export poplgleafmul=1.5
     echo "Default (only relevant to 'core' task): set population delination branch support threshold multiplier for leaf populations to $poplgleafmul"
    fi
    if [ -z "$popbsthresh" ] ; then
     export popbsthresh=80
     echo "Default (only relevant to 'core' task): set population delination branch support threshold to $popbsthresh"
    fi
    if [ -z "$rootingmethod" ] ; then
     export rootingmethod='treebalance'
     echo "Default (only relevant to 'core' task): set reference tree rooting method to $rootingmethod"
    fi
    if [ -z "$hpcremoteptgroot" ] ; then
     echo "Default: all computations will be run locally"
     export hpcremoteptgroot='none'
    fi
}

testmandatoryarg (){
  if [ -z "$2" ]; then
   echo "ERROR: missing argument for option '$1'" 1>&2
   echo "see pantagruel --help for more details" 1>&2
   exit 1
  fi
}

checkexectask (){
  if [ $? != 0 ]; then
    echo "ERROR: Pantagrel pipeline task $1: failed." 1>&2
    exit 1
  else
    echo "Pantagrel pipeline task $1: complete."
  fi
}

promptdate () {
  echo $(date +'[%Y-%m-%d %H:%M:%S]') $1
}

ARGS=`getopt --options "d:r:i:I:f:a:T:A:L:s:t:RV:N:H:cg:hFz" --longoptions "dbname:,rootdir:,initfile:,refresh,iam:,famprefix:,refseq_ass:,refseq_list:,refseq_ass4annot:,refseq_list4annot:,custom_ass:,taxonomy:,pseudocore:,core_seqtype:,pop_lg_thresh:,pop_bs_thresh:,rooting:,reftree:,resume,env_var:,threads:,submit_hpc:,collapse,collapse_param:,genefam_list:,help,FORCE,compress" --name "pantagruel" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  usage ; exit 1
fi

eval set -- "$ARGS"

while true;
do
  case "$1" in
    -h|--help) 
      usagelong
      exit 0;;
    
    -F|--FORCE) 
      export runmode='force'
      shift ;;  
	
    -z|--compress) 
      export compress='on'
      shift ;;
    
    --refresh) 
      export runmode='refreshconfig'
      shift
      export postrefreshargs="${*/--/}"
	  ;;
      
    -d|--dbname) 
      testmandatoryarg "$1" "$2"
      export ptgdbname="$2"
      shift 2;;
    
    -r|--rootdir)
      testmandatoryarg "$1" "$2"
      export ptgroot=$(readlink -f $2)
      shift 2;;
    
    -i|--initfile)
      testmandatoryarg "$1" "$2"
      export initfile=$(readlink -f $2)
      shift 2;;
    
    -I|--iam)
      testmandatoryarg "$1" "$2"
      export myemail="$2"
      echo "set identity to '$myemail'"
      shift 2;;
    
    -f|--famprefix)
      testmandatoryarg "$1" "$2"
      export famprefix="$2"
      echo "set gene family prefix to '$famprefix'"
      shift 2;;

    -a|--custom_ass)
      testmandatoryarg "$1" "$2"
      export customassemb=$(readlink -f $2)
      echo "set custom (raw) genome assembly source folder to '$customassemb'"
      shift 2;;

    -A|--refseq_ass)
      testmandatoryarg "$1" "$2"
      export ncbiass=$(readlink -f $2)
      echo "set NCBI RefSeq(-like) genome assembly source folder to '$ncbiass'"
      shift 2;;

    -L|--refseq_list)
      testmandatoryarg "$1" "$2"
      export listncbiass=$(readlink -f $2)
      echo "set NCBI RefSeq(-like) genome assembly id list to '$listncbiass'"
      shift 2;;

    --refseq_ass4annot)
      testmandatoryarg "$1" "$2"
      export refass=$(readlink -f $2)
      echo "set NCBI RefSeq(-like) genome assembly source folder for reference in user genome annotation to '$refass'"
      shift 2;;

    --refseq_list4annot)
      testmandatoryarg "$1" "$2"
      export listrefass=$(readlink -f $2)
      echo "set NCBI RefSeq(-like) genome assembly id list for reference in user genome annotation to '$listrefass'"
      shift 2;;

    -s|--pseudocore)
      testmandatoryarg "$1" "$2"
      export pseudocoremingenomes="$2"
      echo "set min number of genomes for inclusion in pseudo-core gene set as $pseudocoremingenomes"
      shift 2;;

    -t|--reftree)
      testmandatoryarg "$1" "$2"
      export userreftree=$(readlink -f $2)
      echo "set reference tree as $userreftree"
      shift 2;;

    --core_seqtype)
      testmandatoryarg "$1" "$2"
      export coreseqtype="$2"
      echo "set core genome seuquence type to $coreseqtype"
      case "${coreseqtype}" in
        prot|cds)
         echo "set core genome sequence type to $coreseqtype" ;;
        *)
          echo "ERROR: incorrect core sequence type was specified: ${coreseqtype} (must be 'prot' or 'cds'); exit now"
          exit 1 ;;
      esac
      shift 2;;
    
    --pop_lg_thresh)
      testmandatoryarg "$1" "$2"
      export poplgthresh="$2"
      echo "set population delination branch length threshold to $poplgthresh"
      shift 2;;
    
    --pop_bs_thresh)
      testmandatoryarg "$1" "$2"
      export popbsthresh="$2"
      echo "set population delination branch support threshold to $popbsthresh"
      shift 2;;
      
    --rooting)
      testmandatoryarg "$1" "$2"
      export rootingmethod="$2"
      echo "set reference tree rooting method to $rootingmethod"
      shift 2;;

    -R|--resume)
      export resumetask=true
      echo "will try and resume computation of task where it was last stopped"
      shift ;;

    -V|--env_var)
      testmandatoryarg "$1" "$2"
      export extravars="${2}"
      echo "will add the following environment variable definition to the configuration file: ${extravars}"
      shift 2;;
	  
	-N|--threads)
      testmandatoryarg "$1" "$2"
	  export ptgthreads=${2}
      echo "will try and execute processes in parallel with the following number of threads: ${ptgthreads}"
      shift 2;;
	  
    -T|--taxonomy)
      testmandatoryarg "$1" "$2"
      export ncbitax=$(readlink -f $2)
      echo "set NCBI Taxonomy source folder to '$ncbitax'"
      shift 2;;

    -H|--submit_hpc)
      testmandatoryarg "$1" "$2"
      export hpcremoteptgroot="$2"
      echo "set address of database root folder on remote HPC cluster server"
      shift 2;;

    -c|--collapse)
      export chaintype="collapsed"
      echo "gene tree clade collapsing enabled"
      shift ;;

    --collapse_param)
      testmandatoryarg "$1" "$2"
      export collapseCladeParams="$2"
      echo "set parameters for rake clade collapsing in gene trees"
      shift 2;;

    -g|--genefam_list)
      testmandatoryarg "$1" "$2"
      export genefamlist=$(readlink -f $2)
      echo "set resticted list for computation of gene trees"
      shift 2;;

    --)
      shift
      break;;
    
    *)
  esac
done

tasks=""
while [[ ! -z "${@}" ]] ; do
 case "$1" in
 init)
          tasks="init"
          break ;;
 all)                                     
          tasks="0 1 2 3 4 5 6 7 8 9"
          break;;
 0|00|fetch|fetch_data)
          tasks="${tasks} 0" ;;
 1|01|homologous|homologous_seq_families )
          tasks="${tasks} 1" ;;
 2|02|align|align_homologous_seq)
          tasks="${tasks} 2" ;;
 3|03|sqldb|create_sqlite_db    )
          tasks="${tasks} 3" ;;
 4|04|functional|functional_annotations    )
          tasks="${tasks} 4" ;;
 5|05|core|core_genome_ref_tree)
          tasks="${tasks} 5" ;;
 6|06|genetrees|gene_trees)
          tasks="${tasks} 6" ;;
 7|07|reconciliations)
          tasks="${tasks} 7" ;;
 8|08|specific|clade_specific_genes)
          tasks="${tasks} 8" ;;
 9|09|coevolution)
          tasks="${tasks} 9" ;;
  esac
  shift
done

# sort tasks
tasks=$(echo "${tasks}" | tr ' ' '\n' | sort -u | xargs)
echo "# will run tasks: $tasks"

for task in ${tasks} ; do
  if [[ "$task" == 'init' ]] ; then
    ## init task
    promptdate "Pantagrel pipeline task $task: initiate pangenome database."
    if [ "$runmode" == 'refreshconfig' ] ; then
      if [ -z "$initfile" ] ; then
        echo "Error: cannot use --refresh option without specifying a source config file with -i option"
        exit 1
      fi
      echo "Will refresh the database settings from the configuration file '${initfile}'"
      source ${initfile}
      echo "re-run previous command line:"
      unset initfile
      export runmode='wasrefresh'
      if [ ! -z "${postrefreshargs}" ] ; then
	    echo "(appending the new arguments: ${postrefreshargs})"
        ptginitcmd=${ptginitcmd/init/${postrefreshargs}}
      fi
      echo "# ${ptginitcmd}"
      eval "${ptginitcmd}"
    else
      # check presence of mandatory arguments
      if [ -z "$ptgdbname" ] ; then
       echo -e "Error: Must specify database name\n"
        usage
       exit 1
      fi
      if [ -z "$ptgroot" ] ; then
       echo -e "Error: Must specify root directory location\n"
       usage
       exit 1
      fi
      if [[ -z "$ncbiass" && -z "$listncbiass" && -z "$customassemb" ]] ; then
       echo -e "Error: Must specify at least one folder of input assemblies with options '-A', '-L' or '-a', or any combination of them.\n"
       usage
       exit 1
      fi
      ptgenvsetdefaults
      export ptginitcmd="pantagruel $ptgcmdargs"
      ${ptgscripts}/pipeline/pantagruel_pipeline_init.sh
      checkexectask "$task"
    fi
  else
    ## runtime task
    if [ -z ${initfile} ] ; then
      echo "Error: a Pantagruel configuration file must be provided through option -i ; exit now."
      usage
      exit 1
    fi
    case "$task" in
    0)
     promptdate "Pantagrel pipeline task $task: fetch public genome data from NCBI sequence databases and annotate private genomes."
     ${ptgscripts}/pipeline/pantagruel_pipeline_00_fetch_data.sh ${initfile}
     checkexectask "$task" ;;
    1)
     promptdate "Pantagrel pipeline task $task: classify protein sequences into homologous families."
     ${ptgscripts}/pipeline/pantagruel_pipeline_01_homologous_seq_families.sh ${initfile}
     checkexectask "$task"  ;;
    2)
     promptdate "Pantagrel pipeline task $task: align homologous protein sequences and translate alignemnts into coding sequences."
     ${ptgscripts}/pipeline/pantagruel_pipeline_02_align_homologous_seq.sh ${initfile}
     checkexectask "$task"  ;;
    3)
     promptdate "Pantagrel pipeline task $task: initiate SQL database and load genomic object relationships."
     ${ptgscripts}/pipeline/pantagruel_pipeline_03_create_sqlite_db.sh ${initfile}
     checkexectask "$task"  ;;
    4)
     promptdate "Pantagrel pipeline task $task: use InterProScan to functionally annotate proteins in the database."
     ${ptgscripts}/pipeline/pantagruel_pipeline_04_functional_annotation.sh ${initfile}
     checkexectask "$task"  ;;
    5)
     promptdate "Pantagrel pipeline task $task: select core-genome markers and compute reference tree."
     ${ptgscripts}/pipeline/pantagruel_pipeline_05_core_genome_ref_tree.sh ${initfile}
     checkexectask "$task"  ;;
    6)
     promptdate "Pantagrel pipeline task $task: compute gene trees."
     ${ptgscripts}/pipeline/pantagruel_pipeline_06_gene_trees.sh ${initfile}
     checkexectask "$task"  ;;
    7)
     promptdate "Pantagrel pipeline task $task: compute species tree/gene tree reconciliations."
     ${ptgscripts}/pipeline/pantagruel_pipeline_07_reconciliations.sh ${initfile}
     checkexectask "$task"  ;;
    8)
     promptdate "Pantagrel pipeline task $task: classify genes into orthologous groups (OGs) and search clade-specific OGs."
     ${ptgscripts}/pipeline/pantagruel_pipeline_08_clade_specific_genes.sh ${initfile}
     checkexectask "$task"  ;;
    9)
     promptdate "Pantagrel pipeline task $task: evaluate gene co-evolution and build gene association network."
     ${ptgscripts}/pipeline/pantagruel_pipeline_09_coevolution.sh ${initfile}
     checkexectask "$task"  ;;
     esac
  fi
done
