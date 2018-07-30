#!/bin/bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################

# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 30 July 2018

while getopts u:d:p:f: option
do
case "${option}"
in
u) USER=${OPTARG};;
d) DATE=${OPTARG};;
p) PRODUCT=${OPTARG};;
f) FORMAT=$OPTARG;;
esac
done


export raproot=$1
envsourcescript=${raproot}/environ_pantagruel.sh
source $envsourcescript

shift
tasks=($@)

for task in tasks ; do
 case "$task" in
 
 "init"             ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_init.sh ;;
 "00|fetch_data"    ${ptgscripts}/itemized_pipeline/pantagruel_pipeline_00_fetch_data.sh ;;
 
esac
