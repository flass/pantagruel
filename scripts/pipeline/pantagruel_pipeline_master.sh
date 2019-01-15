#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################
# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 01 October 2018.

usage (){
  echo "Usage: pantagruel -d db_name [-r root_dir] [other options] TASK1 [task-specific options]"
  echo ""
  echo "  _mandatory options_"
  echo "    -d|--dbname       database name"
  echo "    -r|--rootdir      root directory where to create the database; defaults to current folder"
  echo ""
  echo "  _facultative options_"
  echo "    -p|--ptgrepo      location of pantagruel software head folder; defaults to directory where this script is located"
  echo "    -i|--iam          database creator identity (e-mail address is preferred)"
  echo "    -f|--famprefix    alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; defaults to 'PANTAG'"
  echo "                       the chosen prefix will be appended with a 'P' for protein families and a 'C' for CDS families."
  echo "    -T|--taxonomy      path to folder of taxonomy database flat files; defaults to \$rootdir/NCBI/Taxonomy"
  echo "                       if this is not containing the expected file, triggers downloading the daily dump from NCBI Taxonomy at task 00"
  echo "    -A|--refseq_ass  path to folder of source genome assembly flat files formated like NCBI Assembly RefSeq whole directories;"
  echo "                       these can be obtained by searching https://www.ncbi.nlm.nih.gov/assembly and downloadingresults with options:"
  echo "                         Source Database = 'RefSeq' and File type = 'All file types (including assembly-structure directory)'."
  echo "                       defaults to \$rootdir/NCBI/Assembly"
  echo "    -a|--custom_ass  path to folder of user-provided genomes, containing:"
  echo "                      _mandatory_ "
  echo "                       - a 'contigs/' folder where are stored all source genome assembly FASTA files"
  echo "                           OR"
  echo "                       - a 'prokka_annotation/' folder where are stored all files resulting from Prokka annotation"
  echo "                      _optionally_ "
  echo "                       - a 'strain_infos.txt' file"
  echo "                       (unnanotated contig fasta files); defaults to \$rootdir/user_genomes"
  echo "    -s|--pseudocore  integer number, the minimum number of genomes in which a gene family should be present"
  echo "                       to be included in the pseudo-core genome gene set (otherwise has to be set interactively"
  echo "                       when running task 'core')."
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
  echo "    -c|--collapse      enable collapsing the rake clades in the gene trees (strongly recomended in datasets of size > 50 genomes)."
  echo "    -C|--collapse_par  [only for 'genetrees' task] specify parameters for collapsing the rake clades in the gene trees."
  echo "                       A single-quoted, semicolon-delimited string containing variable definitions must be provided."
  echo "                       Default is equivalent to providing the following string:"
  echo "                          'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'"
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
  if [ -z "$2" ]; then echo "missing argument for option '$1'" usage ; exit 1 ; fi
}
checkexec (){
  if [ $? != 0 ]; then echo "ERROR: last pipeline step failed" ; exit 1 ; fi
}


ARGS=`getopt --options "d:r:p:i:f:a:T:A:s:H:cC:h" --longoptions "dbname:,rootdir:,ptgrepo:,iam:,famprefix:,refseq_ass:,custom_ass:,taxonomy:,pseudocore:,submit_hpc:,collapse,collapse_par:,help" --name "pantagruel" -- "$@"`

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

    -T|--taxonomy)
      testmandatoryarg "$1" "$2"
      export pseudocoremingenomes="$2"
      echo "set min number of genomes for inclusion in pseudo-core gene set"
      shift 2;;

    -s|--pseudocore)
      testmandatoryarg "$1" "$2"
      export ncbitax="$2"
      echo "set NCBI Taxonomy source folder to '$ncbitax'"
      shift 2;;

    -H|--submit_hpc
      testmandatoryarg "$1" "$2"
      export hpcremoteptgroot="$2"
      echo "set address of database root folder on remote HPC cluster server"
      shift 2;;

    -c|--collapse
      export chaintype="collapsed"
      echo "gene tree clade collapsing enabled"
      shift ;;

    -C|--collapse_par
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
 echo "Must specify database name"
 usage
 exit 1
fi
if [ -z "$ptgroot" ] ; then
 echo "Must specify root directory location"
 usage
 exit 1
fi


# Default values:
if [ -z "$ptgrepo" ] ; then
  export ptgrepo=$(dirname $0)
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
  if [ ! -d "$customassemb/contigs" ] ; then
   echo "custom assembly folder must contain a subfolder 'contigs/'"
   usage
   exit 1
  fi
  echo "Default: set custom (raw) genome assembly source folder to '$customassemb'"
fi
if [ -z "$ncbiass" ] ; then
  export ncbiass=${ptgroot}/NCBI/Assembly
  echo "Default: set NCBI RefSeq(-like) genome assembly source folder to '$ncbiass'"
fi
if [ -z "$ncbitax" ] ; then
 export ncbitax=${ptgroot}/NCBI/Taxonomy
 echo "Default: set NCBI Taxonomy source folder to '$ncbitax'"
fi

if [ -z "$hpcremoteptgroot" ] ; then
 export hpcremoteptgroot='none'
 echo "Default: all computations will be run locally (instead of on a potential HPC service)"
fi

if [ -z "$chaintype" ] ; then
 export chaintype="fullgenetree"
 echo "Default: gene tree clade collapsing disabled"
fi

if [ -z "$collapseCladeParams" ] ; then
  export collapseCladeParams='default'
  if [ ! -z "$chaintype" ] ; then 
    echo "Default: gene tree clade collapsing will use these parameters: 'cladesupp=70 ; subcladesupp=35 ; criterion=bs ; withinfun=median'"
  fi
fi

echo "will create/use Pantagruel database '$ptgdbname', set in root folder: '$ptgroot'"
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
          tasks="${tasks} 3" ;;
 5|05|core|core_genome_ref_tree)
          tasks="${tasks} 4" ;;
 6|06|genetrees|gene_trees)
          tasks="${tasks} 5" ;;
 7|07|reconciliations)
          tasks="${tasks} 6" ;;
 8|08|specific|clade_specific_genes)
          tasks="${tasks} 8" ;;
 9|09|coevolution)
          tasks="${tasks} 7" ;;
  esac
  shift
done

# sort tasks
tasks=$(echo "${tasks}" | tr ' ' '\n' | sort -u | xargs)
echo $tasks

for task in "$tasks" ; do
  if [[ "$task" == 'init' ]] ; then
    echo "Pantagrel pipeline step $task: initiate pangenome database."
    ${ptgscripts}/pipeline/pantagruel_pipeline_init.sh ${ptgdbname} ${ptgroot} ${ptgrepo} ${myemail} ${famprefix} ${ncbiass} ${ncbitax} ${customassemb} ${hpcremoteptgroot} ${initfile}
    checkexec
  else
   export ptgdb=${ptgroot}/${ptgdbname}
   envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
   source $envsourcescript
   case "$task" in
   0)
    echo "Pantagrel pipeline step $task: fetch public genome data from NCBI sequence databases and annotate private genomes."
    ${ptgscripts}/pipeline/pantagruel_pipeline_00_fetch_data.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   1)
    echo "Pantagrel pipeline step $task: classify protein sequences into homologous families."
    ${ptgscripts}/pipeline/pantagruel_pipeline_01_homologous_seq_families.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   2)
    echo "Pantagrel pipeline step $task: align homologous protein sequences and translate alignemnts into coding sequences."
    ${ptgscripts}/pipeline/pantagruel_pipeline_02_align_homologous_seq.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   3)
    echo "Pantagrel pipeline step $task: initiate SQL database and load genomic object relationships."
    ${ptgscripts}/pipeline/pantagruel_pipeline_03_create_sqlite_db.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   4)
    echo "Pantagrel pipeline step $task: use InterProScan to functionally annotate proteins in the database."
    ${ptgscripts}/pantagruel_pipeline_04_functional_annotation.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   5)
    echo "Pantagrel pipeline step $task: select core-genome markers and compute reference tree."
    if [[ -z ${pseudocoremingenomes} || "${pseudocoremingenomes}"=='REPLACEpseudocoremingenomes' ]] ; then
      echo "'pseudocoremingenomes' variable is not set; will run INTERACTIVELY $ptgscripts/choose_min_genome_occurrence_pseudocore_genes.sh to choose a sensible value."
      echo ""
      ${ptgscripts}/choose_min_genome_occurrence_pseudocore_genes.sh
      echo ""
    fi
    ${ptgscripts}/pipeline/pantagruel_pipeline_05_core_genome_ref_tree.sh ${ptgdbname} ${ptgroot} ${pseudocoremingenomes}
    checkexec  ;;
   6)
    echo "Pantagrel pipeline step $task: compute gene trees."
    ${ptgscripts}/pipeline/pantagruel_pipeline_06_gene_trees.sh ${ptgdbname} ${ptgroot} ${collapseCladeParams} ${hpcremoteptgroot}
    checkexec  ;;
   7)
    echo "Pantagrel pipeline step $task: compute species tree/gene tree reconciliations."
    ${ptgscripts}/pipeline/pantagruel_pipeline_07_reconciliations.sh ${ptgdbname} ${ptgroot} ${hpcremoteptgroot}
    checkexec  ;;
   8)
    echo "Pantagrel pipeline step $task: classify genes into orthologous groups (OGs) and search clade-specific OGs."
    ${ptgscripts}/pipeline/pantagruel_pipeline_08_clade_specific_genes.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
   9)
    echo "Pantagrel pipeline step $task: evaluate gene co-evolution and build gene association network."
    ${ptgscripts}/pipeline/pantagruel_pipeline_09_coevolution.sh ${ptgdbname} ${ptgroot}
    checkexec  ;;
    esac
  fi
done
