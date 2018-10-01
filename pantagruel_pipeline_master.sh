#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################
# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 01 October 2018.

usage (){
  echo "Usage: pantagruel_pipeline_master.sh -d db_name [-r root_dir] [other options] TASK1 [TASK2 [...]]"
  echo "    -d|-dbname     database name"
  echo "    -r|-rootdir    root directory where to create the database; defaults to current folder"
  echo "    -p|-ptgrepo    location of pantagruel software head folder; default to directory where this script is located"
  echo "    -i|-iam        database creator identity (e-mail address is preferred)"
  echo "    -f|-famprefix  alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters;"
  echo "                   the chosen prefix will be appended with a 'P' for protein families and a 'C' for CDS families."
  echo "TASKs are to be picked among the following (equivalent digit/number/keywords are separated by a '|'):"
  echo "  all"
  echo "       perform all pipeline tasks, from database initiation to co-evolution inference"
  echo "  init"
  echo "       database initiation: define shared environment variables"
  echo "  0|00|fetch|fetch_data"  
  echo "       fetch public genome data from NCBI sequence databases and annotate private genomes"
  echo "  1|01|homologous|homologous_seq_families"
  echo "       classify protein sequences into homologous families"
  echo "  2|02|align|align_homologous_seq"  
  echo "       align homologous protein sequences and translate alignemnts into coding sequences"
  echo "  3|03|sqldb|create_sqlite_db"    
  echo "       initiate SQL database and load genomic object relationships"  
  echo "  4|04|core|core_genome_ref_tree"
  echo "       select core-genome markers and compute reference tree"
  echo "  5|05|genetrees|gene_trees"
  echo "       compute gene tree"
  echo "  6|06|reconciliations"
  echo "       compute species tree/gene tree reconciliations"
  echo "  7|07|coevolution"
  echo "       quantify gene co-evolution and build gene association network"
  echo "  8|08|specific|clade_specific_genes"
  echo "       classify genes into orthologous groups (OGs) and search clade-specific OGs"
}
testmandatoryarg (){
  if [ -z "$2" ]; then echo "missing argument for option '$1'" usage ; exit 1 ; fi
}
checkexec (){
  if [ $? != 0 ]; then echo "ERROR: last pipeline step failed" ; exit 1 ; fi
}


# Default values:
export ptgroot=${PWD}
export ptgrepo=$(dirname $0)
export myemail="undisclosed"
export famprefix="PANTAG"

ARGS=`getopt --options "d:r:p:i:f:" --longoptions "dbname:,rootdir:,ptgrepo:,iam:,famprefix:" --name "pantagruel_pipelin_master.sh" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  usage ; exit 1
fi

eval set -- "$ARGS"

while true;
do
  case "$1" in
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
      shift 2;;
    
    -i|--iam)
      testmandatoryarg "$1" "$2"
      export myemail="$2"
      shift 2;;
    
    -f|--famprefix)
      testmandatoryarg "$1" "$2"
      export famprefix="$2"
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

echo "will create/use Pantagruel database '$ptgdbname' in local root folder: '$PWD'"


## optional parameters
export assembler=""
export seqcentre=""


tasks=""
for arg in $@ ; do
 case "$task" in
 "all")                                        tasks="init 0 1 2 3 4 5 6 7 8" ; break;;
 "init")                                       tasks="init ${tasks}" ;;
 "0|00|fetch|fetch_data"  )                    tasks="${tasks} 0" ;;
 "1|01|homologous|homologous_seq_families" )   tasks="${tasks} 1" ;;
 "2|02|align|align_homologous_seq"    )        tasks="${tasks} 2" ;;
 "3|03|sqldb|create_sqlite_db"    )            tasks="${tasks} 3" ;;
 "4|04|core|core_genome_ref_tree")             tasks="${tasks} 4" ;;
 "5|05|genetrees|gene_trees")                  tasks="${tasks} 5" ;;
 "6|06|reconciliations")                       tasks="${tasks} 6" ;;
 "7|07|coevolution")                           tasks="${tasks} 7" ;;
 "8|08|specific|clade_specific_genes")         tasks="${tasks} 8" ;;
  esac
done

if [[ ${tasks[0]} != "all" && ${tasks[0]} != "init" ]] ; then
  export ptgdb=${ptgroot}/${ptgdbname}
  envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
  source $envsourcescript
fi

for task in "$tasks" ; do
 case "$task" in
 "init")
    echo "Pantagrel pipeline step $task: initiate pangenome database."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_init.sh
    checkexec ;;
 "0")
    echo "Pantagrel pipeline step $task: fetch public genome data from NCBI sequence databases and annotate private genomes."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_00_fetch_data.sh
    checkexec  ;;
 "1")
    echo "Pantagrel pipeline step $task: classify protein sequences into homologous families."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_01_homologous_seq_families.sh
    checkexec  ;;
 "2")
    echo "Pantagrel pipeline step $task: align homologous protein sequences and translate alignemnts into coding sequences."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_02_align_homologous_seq.sh
    checkexec  ;;
 "3")
    echo "Pantagrel pipeline step $task: initiate SQL database and load genomic object relationships."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_03_create_sqlite_db.sh
    checkexec  ;;
 "4")
    echo "Pantagrel pipeline step $task: select core-genome markers and compute reference tree."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_04_core_genome_ref_tree.sh
    checkexec  ;;
 "5")
    echo "Pantagrel pipeline step $task: compute gene trees."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_05_gene_trees.sh
    checkexec  ;;
 "6")
    echo "Pantagrel pipeline step $task: compute species tree/gene tree reconciliations."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_06_reconciliations.sh
    checkexec  ;;
 "7")
    echo "Pantagrel pipeline step $task: evaluate gene co-evolution and build gene association network."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_07_coevolution.sh
    checkexec  ;;
 "8")
    echo "Pantagrel pipeline step $task: classify genes into orthologous groups (OGs) and search clade-specific OGs."
    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_08_clade_specific_genes.sh
    checkexec  ;;
  esac
done
