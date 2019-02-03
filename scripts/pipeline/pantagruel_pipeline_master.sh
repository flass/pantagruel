#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################
# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 01 October 2018.

locexec=$(readlink -e $0)
defptgrepo=${locexec%%/scripts/pipeline/*}

usage (){
  echo "Usage: pantagruel -d db_name [-r root_dir] [other options] TASK1 [task-specific options]"
  echo ""
  echo "  _mandatory options_"
  echo "    -d|--dbname       database name"
  echo ""
  echo "    -r|--rootdir      root directory where to create the database; defaults to current folder"
  echo ""
  echo "  _facultative options_"
  echo "    -p|--ptgrepo      location of pantagruel software head folder; defaults to ${defptgrepo}"
  echo ""
  echo "    -i|--iam          database creator identity (e-mail address is preferred)"
  echo ""
  echo "    -f|--famprefix    alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; defaults to 'PANTAG'"
  echo "                       the chosen prefix will be appended with a 'P' for protein families and a 'C' for CDS families."
  echo ""
  echo "    -T|--taxonomy      path to folder of taxonomy database flat files; defaults to \$rootdir/NCBI/Taxonomy_YYYY-MM-DD (suffix is today's date)"
  echo "                        if this is not containing the expected file, triggers downloading the daily dump from NCBI Taxonomy at task 00"
  echo ""
  echo "    -A|--refseq_ass  path to folder of source genome assembly flat files formated like NCBI Assembly RefSeq whole directories;"
  echo "                       these can be obtained by searching https://www.ncbi.nlm.nih.gov/assembly and downloadingresults with options:"
  echo "                         Source Database = 'RefSeq' and File type = 'All file types (including assembly-structure directory)'."
  echo "                       defaults to \$rootdir/NCBI/Assembly_YYYY-MM-DD (suffix is today's date). "
  echo "                       A simple archive 'genome_assemblies.tar' (as obtained from the NCBI website)can be placed in that folder."
  echo "                       If user genomes are also provided, these RefSeq assemblies will be used as reference for their annotation."
  echo ""
  echo "    --refseq_ass4annot idem, but WILL NOT be used in the study, only as a reference to annotate user genomes (defaults to vaule of -A option)"
  echo ""
  echo "    -a|--custom_ass  path to folder of user-provided genomes (defaults to \$rootdir/user_genomes), containing:"
  echo "                      _mandatory_ "
  echo "                       - a 'contigs/' folder, where are stored multi-FASTA files of genome assemblies (one file per genome,"
  echo "                          with extension '.fa', '.fasta' or '.fas' ...). Fasta file names will be truncated by removing"
  echo "                          the '.fa' string and everything occuring after) and will be retained as the assembly_id (beware "
  echo "                          of names redundant with RefSeq assemblies)."
  echo "                       - a 'strain_infos.txt' file describing the organism, with columns headed:"
  echo "                           'assembly_id'; 'genus'; 'species'; 'strain'; 'taxid'; 'locus_tag_prefix'"
  echo "                         'assembly_id' must match the name of a contig file (e.g. 'assembId.fasta')"
  echo "                         'locus_tag_prefix' must match the prefix of ids given to CDS, proteins and genome regions (contigs)"
  echo "                         in potentially provided annotation files (see below)."
  echo "                      _optional_ "
  echo "                       - an 'annotation/' folder, where are stored annotation files: "
  echo "                         - one mandatory in GFF 3.0 file format (with a '.gff' extension);"
  echo "                          and optionally, the following files (with consistent ids!!):"
  echo "                         - one in GenBank flat file format (with a '.gbk' extension);"
  echo "                         - one in Fasta format containing CDS sequences (with a '.ffn' extension)."
  echo "                         - one in Fasta format containing matching protein sequences (with a '.faa' extension)."
  echo "                         These four files are produced when using Prokka for annotation; if at least one of the .gbk, .ffn or .faa"
  echo "                         are missing, all three will be derived from the .gff source. Each genome annotation file set must be stored"
  echo "                         in a separate folder, which name must match a contig file (e.g. 'assembId/' for 'assembId.fasta')."
  echo "                      NOTE: to ensure proper parsing, it is strongly advised that any provided annotation was generated with Prokka"
  echo "                      NOTE: to ensure uniform annotation of the dataset, it is advised to let Pantagruel annotate the contigs (calling Prokka)"
  echo ""
  echo "    -s|--pseudocore  integer or string, the minimum number of genomes in which a gene family should be present to be included in"
  echo "                       the pseudo-core genome, i.e. the gene set which alignments will be concatenated for reference tree search."
  echo "                       Only relevant when running task 'core'; a non-integer value will trigger an INTERACTIVE prompt for search of an optimal value."
  echo "                       Defaults to the total number of genomes (work with a strict core genome set)."
  echo ""
  echo "    -t|--reftree     specify a reference tree for reconciliation and clade-specific gene analyses;"
  echo "                       over-rides the computation of tree from the concatenate of (pseudo-)core genome gene during taske 'core'."
  echo ""
  echo "    --core_seqtype   {cds|prot} define the type of sequence that will be used to compute the (pseudo-)core genome tree (default to 'cds')"
  echo ""
  echo "    --pop_lg_thresh  definee the threshold of branch length for delinating populations in the reference tree "
  echo "                       (default: 0.0005 for nucleotide alignemnt-based tree; 0.0001 for protein-based)"
  echo ""
  echo "    --pop_bs_thresh  definee the threshold of branch support for delinating populations in the reference tree (default: 80)"
  echo ""
  echo "    -R|--resume      try and resume the task from previous run that was interupted (for the moment only available for taske 'core')"
  echo ""
  echo "    -H|--submit_hpc  full address (hostname:/folder/location) of a folder on a remote high-performance computating (HPC) cluster server"
  echo "                       This indicate that computationally intensive tasks, including building the gene tree collection"
  echo "                       ('genetrees') and reconciling gene tree with species tree ('reconciliations') will be run"
  echo "                       on a HPC server (only Torque/PBS job submission system is supported so far)."
  echo "                       [support for core genome tree building ('core') remains to be implemented]."
  echo "                       Instead of running the computations, scripts for cluster job submission will be generated automatically."
  echo "                       Data and scripts will be transfered to the specified address (the database folder structure"
  echo "                       will be duplicated there, but only relevant files will be synced). Note that job submission"
  echo "                       scripts will need to be executed manually on the cluster server."
  echo "                       If set at init stage, this option will be maintained for all tasks. However, the remote address"
  echo "                       can be updated when calling a specific task; string 'none' cancels the HPC behaviour."
  echo ""
  echo "    -c|--collapse      enable collapsing the rake clades in the gene trees (strongly recomended in datasets of size > 50 genomes)."
  echo ""
  echo "    -C|--collapse_par  [only for 'genetrees' task] specify parameters for collapsing the rake clades in the gene trees."
  echo "                       A single-quoted, semicolon-delimited string containing variable definitions must be provided."
  echo "                       Default is equivalent to providing the following string:"
  echo "                          'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'"
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
  echo "________________________________________________________________________"
  echo ""
}
testmandatoryarg (){
  if [ -z "$2" ]; then
   echo "ERROR: missing argument for option '$1'" 1>&2
   echo "see pantagruel --help for more details" 1>&2
   exit 1
  fi
}

checkexec (){
  if [ $? != 0 ]; then
    echo "ERROR: $1" 1>&2
    exit 1
  else
    if [ ! -z "$2" ] ; then
      echo "$2"
    fi
  fi
}

promptdate () {
  echo $(date +'[%Y-%m-%d %H:%M:%S]') $1
}

ARGS=`getopt --options "d:r:p:i:f:a:T:A:s:t:RH:cC:h" --longoptions "dbname:,rootdir:,ptgrepo:,iam:,famprefix:,refseq_ass:,refseq_ass4annot:,custom_ass:,taxonomy:,pseudocore:,core_seqtype:,pop_lg_thresh:,pop_bs_thresh:,reftree:,resume,submit_hpc:,collapse,collapse_par:,help" --name "pantagruel" -- "$@"`

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
      usage
      exit 0;;
    
    -d|--dbname) 
      testmandatoryarg "$1" "$2"
      export ptgdbname="$2"
      shift 2;;
    
    -r|--rootdir)
      testmandatoryarg "$1" "$2"
      export ptgroot="$2"
      shift 2;;
    
    -p|--ptgrepo)
      testmandatoryarg "$1" "$2"
      export ptgrepo="$2"
      echo "set Pantagruel software repository to '$ptgrepo'"
      shift 2;;
    
    -i|--iam)
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
      export customassemb="$2"
      echo "set custom (raw) genome assembly source folder to '$customassemb'"
      shift 2;;

    -A|--refseq_ass)
      testmandatoryarg "$1" "$2"
      export ncbiass="$2"
      echo "set NCBI RefSeq(-like) genome assembly source folder to '$ncbiass'"
      shift 2;;

    --refseq_ass4annot)
      testmandatoryarg "$1" "$2"
      export refass="$2"
      echo "set NCBI RefSeq(-like) genome assembly source folder for reference in user genome annotation to '$refass'"
      shift 2;;

    -s|--pseudocore)
      testmandatoryarg "$1" "$2"
      export pseudocoremingenomes="$2"
      echo "set min number of genomes for inclusion in pseudo-core gene set as $pseudocoremingenomes"
      shift 2;;

    -t|--reftree)
      testmandatoryarg "$1" "$2"
      export userreftree="$2"
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

    -R|--resume)
      export resumetask=true
      echo "will try and resume computation of task where it was last stopped"
      shift ;;

    -T|--taxonomy)
      testmandatoryarg "$1" "$2"
      export ncbitax="$2"
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

    -C|--collapse_par)
      testmandatoryarg "$1" "$2"
      export collapseCladeParams="$2"
      echo "set parameters for rake clade collapsing in gene trees"
      shift 2;;

    --)
      shift
      break;;
    
    *)
  esac
done

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

setdefaults (){
    # Default values:
    if [ -z "$ptgrepo" ] ; then
    export ptgrepo=$defptgrepo
    echo "Default: set Pantagruel software repository to '$ptgrepo'"
    fi
    if [ -z "$myemail" ] ; then
    export myemail="undisclosed"
    echo "Default: set identity to '$myemail'"
    fi
    if [ -z "$famprefix" ] ; then
    export famprefix="PANTAG"
    echo "Default: set gene family prefix to '$famprefix'"
    fi
    if [ -z "$customassemb" ] ; then
    export customassemb=${ptgroot}/user_genomes
    echo "Default: set custom (raw) genome assembly source folder to '$customassemb'"
    fi
    if [ -z "$ncbiass" ] ; then
    export ncbiass=${ptgroot}/NCBI/Assembly_$(date +'%Y-%m-%d')
    echo "Default: set NCBI RefSeq(-like) genome assembly source folder to '$ncbiass'"
    fi
    if [ -z "$refass" ] ; then
    export refass=${ncbiass}
    echo "Default: set NCBI RefSeq(-like) genome assembly source folder for reference in user genome annotation to '$refass'"
    fi
    if [ -z "$ncbitax" ] ; then
     export ncbitax=${ptgroot}/NCBI/Taxonomy_$(date +'%Y-%m-%d')
     echo "Default: set NCBI Taxonomy source folder to '$ncbitax'"
    fi
}

setnondefaults (){
	nondefvardecl=${ptgtmp}/nondefvardecl.sh
	rm -f ${ptgtmp}/nondefvardecl.sh
    if [ ! -z "$hpcremoteptgroot" ] ; then
     echo "Overide default: all computations will be run on a distant server (potentially a HPC service) instead of locally"
     echo "export hpcremoteptgroot=${hpcremoteptgroot}" >> ${ptgtmp}/nondefvardecl.sh
    fi
    if [ ! -z "$chaintype" ] ; then
     echo "Overide default: gene tree clade collapsing is enabled"
     echo "export chaintype=${chaintype}" >> ${ptgtmp}/nondefvardecl.sh
    fi
    if [ ! -z "$collapseCladeParams" ] ; then
     if [ ! -z "$chaintype" ] ; then 
      echo "Overide default: gene tree clade collapsing will use these parameters: '$collapseCladeParams'"
     echo "export collapseCladeParams=${collapseCladeParams}" >> ${ptgtmp}/nondefvardecl.sh
     fi
    fi
    #~ if [ ! -z "$pseudocoremingenomes" ] ; then
     #~ echo "Overide default: will establish pseudo-core gene set (gene present in a minimum of $pseudocoremingenomes genome) instead of a strict core-genome"
     #~ echo "export pseudocoremingenomes=${pseudocoremingenomes}" >> ${ptgtmp}/nondefvardecl.sh
    #~ fi
    if [ ! -z "$userreftree" ] ; then
     echo "Overide default (only relevant to 'core' task): a reference tree was provided; will not compute core-genome tree"
     echo "export userreftree=${userreftree}" >> ${ptgtmp}/nondefvardecl.sh
    fi
    if [ ! -z "$coreseqtype" ] ; then
     echo "Overide default (only relevant to 'core' task): core sequence type is set to '$coreseqtype'"
     echo "export coreseqtype=${coreseqtype}" >> ${ptgtmp}/nondefvardecl.sh
    fi
    if [ ! -z "$poplgthresh" ] ; then
     echo "Overide default (only relevant to 'core' task): set population delination branch support threshold to $poplgthresh"
     echo "export poplgthresh=${poplgthresh}" >> ${ptgtmp}/nondefvardecl.sh
    fi
    if [ ! -z "$popbsthresh" ] ; then
     echo "Overide default (only relevant to 'core' task): set population delination branch support threshold to $popbsthresh"
     echo "export popbsthresh=${popbsthresh}" >> ${ptgtmp}/nondefvardecl.sh
    fi
    if [ ! -z "$resumetask" ] ; then
     echo "Overide default (only relevant to 'core' task): will try and resume computation of task where it was last stopped"
     echo "export resumetask=${resumetask}" >> ${ptgtmp}/nondefvardecl.sh
    fi
}

echo -e "# will create/use Pantagruel database '$ptgdbname', set in root folder: '$ptgroot'\n#"
export ptgscripts=${ptgrepo}/scripts

## task-specific parameters
initfile=""

tasks=""
while [[ ! -z "${@}" ]] ; do
 case "$1" in
 init)
          tasks="init"
          if [ ! -z "$2" ] ; then
             if [[ "$2"=='help' ]] ; then
               echo "Usage: pantagruel [general options] init [init_file]"
               exit 0
             else
               initfile="$2"
             fi
          fi
          break ;;
 all)                                     
          tasks="0 1 2 3 4 5 6 7 8"
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
echo $tasks

for task in "$tasks" ; do
  if [[ "$task" == 'init' ]] ; then
    promptdate "Pantagrel pipeline step $task: initiate pangenome database."
    setdefaults
    ${ptgscripts}/pipeline/pantagruel_pipeline_init.sh ${ptgdbname} ${ptgroot} ${ptgrepo} ${myemail} ${famprefix} ${ncbiass} ${ncbitax} ${customassemb} ${refass} ${initfile}
    checkexec
  else
   case "$task" in
   0)
    promptdate "Pantagrel pipeline step $task: fetch public genome data from NCBI sequence databases and annotate private genomes."
    ${ptgscripts}/pipeline/pantagruel_pipeline_00_fetch_data.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   1)
    promptdate "Pantagrel pipeline step $task: classify protein sequences into homologous families."
    ${ptgscripts}/pipeline/pantagruel_pipeline_01_homologous_seq_families.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   2)
    promptdate "Pantagrel pipeline step $task: align homologous protein sequences and translate alignemnts into coding sequences."
    ${ptgscripts}/pipeline/pantagruel_pipeline_02_align_homologous_seq.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   3)
    promptdate "Pantagrel pipeline step $task: initiate SQL database and load genomic object relationships."
    ${ptgscripts}/pipeline/pantagruel_pipeline_03_create_sqlite_db.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   4)
    promptdate "Pantagrel pipeline step $task: use InterProScan to functionally annotate proteins in the database."
    ${ptgscripts}/pipeline/pantagruel_pipeline_04_functional_annotation.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   5)
    promptdate "Pantagrel pipeline step $task: select core-genome markers and compute reference tree."
    if [[ ! -z "${pseudocoremingenomes}" ]] ; then
	  case "${pseudocoremingenomes}" in
        ''|*[!0-9]*)
          echo "'pseudocoremingenomes' variable is not set to a correct integer value: '${pseudocoremingenomes}' ; unset this variable"
          unset pseudocoremingenomes
          echo "will run INTERACTIVELY '$ptgscripts/choose_min_genome_occurrence_pseudocore_genes.sh' to choose a sensible value."
          echo ""
          ${ptgscripts}/choose_min_genome_occurrence_pseudocore_genes.sh ${ptgdbname} ${ptgroot} ;;
        *)
          echo "'pseudocoremingenomes' variable is set to ${pseudocoremingenomes}"
          echo "will run non-interactively '$ptgscripts/choose_min_genome_occurrence_pseudocore_genes.sh' to record the gene family set."
          echo ""
          ${ptgscripts}/choose_min_genome_occurrence_pseudocore_genes.sh ${ptgdbname} ${ptgroot} ${pseudocoremingenomes} ;;
      esac
    else
       echo "'pseudocoremingenomes' variable is not set"
       echo "will rely on strict core genome definition"
       echo "will run non-interactively $ptgscripts/choose_min_genome_occurrence_pseudocore_genes.sh to record the gene family set."
       echo ""
       ${ptgscripts}/choose_min_genome_occurrence_pseudocore_genes.sh ${ptgdbname} ${ptgroot} 'REPLACEpseudocoremingenomes'
    fi
    setnondefaults
    ${ptgscripts}/pipeline/pantagruel_pipeline_05_core_genome_ref_tree.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   6)
    promptdate "Pantagrel pipeline step $task: compute gene trees."
    setnondefaults
    ${ptgscripts}/pipeline/pantagruel_pipeline_06_gene_trees.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   7)
    promptdate "Pantagrel pipeline step $task: compute species tree/gene tree reconciliations."
    setnondefaults
    ${ptgscripts}/pipeline/pantagruel_pipeline_07_reconciliations.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   8)
    promptdate "Pantagrel pipeline step $task: classify genes into orthologous groups (OGs) and search clade-specific OGs."
    ${ptgscripts}/pipeline/pantagruel_pipeline_08_clade_specific_genes.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   9)
    promptdate "Pantagrel pipeline step $task: evaluate gene co-evolution and build gene association network."
    ${ptgscripts}/pipeline/pantagruel_pipeline_09_coevolution.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
    esac
  fi
  rm -f ${ptgtmp}/nondefvardecl.sh
done
