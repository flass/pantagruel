#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################
# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

usage (){
  echo "Usage: pantagruel_pipeline_master.sh [options] [db_name] [root_dir]"
  echo "    -d|-dbname     database name; overrides non option argument $1"
  echo "    -r|-rootdir    root directory where to create the database; overrides non option argument $2; defaults to current folder"
  echo "    -p|-ptgrepo    location of pantagruel software head folder; default to directory where this script is located"
  echo "    -i|-iam        database creator identity (e-mail address is preferred)"
  echo "    -f|-famprefix  alphanumerical prefix (no number first) of the names for homologous protein/gene family clusters; will be appended with a 'P' for proteins and a 'C' for CDSs."
}
testmandatoryarg (){
  if [ -z "$2" ]; then echo "missing argument for option '$1'" usage ; exit 1 ; fi
}


today=$(date +%d-%m-%Y)
export ptgrepo=$(dirname $0)
export myemail="me.myself@respectable-institu.ti.on"
export ptgroot=${PWD}
export famprefix="PANTAG"
export ptgdbname="PantagruelDB_${today}"
export famprefix="ABCDE"

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
  if [ -n "$1" ] ; then
    export ptgdbname="$1"
    shift
  else
    echo "Must specify database name"
    usage
    exit 1
fi
if [ -z "$ptgroot" ] ; then
  if [ -n "$1" ] ; then
    export ptgroot="$1"
    shift
  else
    echo "Must specify root directory in which to create database"
    usage
    exit 1
fi


## optional parameters
export assembler=""
export seqcentre=""


tasks=""
for arg in $@ ; do
 case "$task" in
 "init"                                       tasks="init ${tasks}" ;;
 "all"                                        tasks="init 0 1 2 3 4 5 6 7 8" ;;
 "0|00|fetch|fetch_data"                      tasks="${tasks} 0" ;;
 "1|01|homologous|homologous_seq_families"    tasks="${tasks} 1" ;;
 "2|02|align|align_homologous_seq"            tasks="${tasks} 2" ;;
 "3|03|sqldb|create_sqlite_db"                tasks="${tasks} 3" ;;
 "4|04|core|core_genome_ref_tree"             tasks="${tasks} 4" ;;
 "5|05|genetrees|gene_trees"                  tasks="${tasks} 5" ;;
 "6|06|reconciliations"                       tasks="${tasks} 6" ;;
 "7|07|coevolution"                           tasks="${tasks} 7" ;;
 "8|08|specific|clade_specific_genes"         tasks="${tasks} 8" ;;

  esac
done

if [[ ${tasks[0]} != "all" && ${tasks[0]} != "init" ]] ; then
  export ptgdb=${ptgroot}/${ptgdbname}
  envsourcescript=${ptgdb}/environ_pantagruel_${ptgdbname}.sh
  source $envsourcescript
fi

for task in "$tasks" ; do
 case "$task" in
 "init"                                       ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_init.sh ;;
 "0"                      ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_00_fetch_data.sh ;;
 "1"    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_01_homologous_seq_families.sh ;;
 "2"            ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_02_align_homologous_seq.sh ;;
 "3"                ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_03_create_sqlite_db.sh ;;
 "4"             ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_04_core_genome_ref_tree.sh ;;
 "5"                  ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_05_gene_trees.sh ;;
 "6"                       ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_06_reconciliations.sh ;;
 "7"                           ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_07_coevolution.sh ;;
 "8"         ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_08_clade_specific_genes.sh ;;

  esac
done
