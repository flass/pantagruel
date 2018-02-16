#!/bin/bash

assembly_folder=$1
DEST=$2

assemblies=$(ls ${assembly_folder}/GC[AF]_* -d)
faas=$(echo $assemblies | sed -e 's#\([^ ]\+\)/\([^/ ]\+\)#\1/\2/\2_protein\.faa\.gz#g')
for faa in $faas; do zcat $faa ; done > ${DEST}
