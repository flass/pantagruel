#!/bin/bash

tasklist=$1
outdir=$2
nbthreads=2
if [ ! -z $3 ] ; then
  nbthreads=$3
fi
overwrite=0
if [ $4 == 'overwrite' ] ; then
  overwrite=1
fi

clustalo=/usr/bin/clustalo

for task in `cat $tasklist` ; do
  bn=`basename $task`
  rad=${bn%.fasta}
  if [[ ! -e $outdir/$rad.aln || $overwrite -eq 1 ]] ; then
    echo "task: $task"
    date +"%d/%m/%Y %H:%M:%S"
    clustalo --threads=$nbthreads -i $task -o $outdir/$rad.aln
    date +"%d/%m/%Y %H:%M:%S"
    echo ""
    echo "- - - - -"
    echo ""
  fi
done

