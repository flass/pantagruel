#!/bin/bash

assembly_list=$1
DEST=$2

faas=$(cat $assembly_list | sed -e 's#\([^ ]\+\)/\([^/ ]\+\)#\1/\2/\2_protein\.faa\.gz#g')
for faa in $faas; do zcat $faa ; done > ${DEST}
