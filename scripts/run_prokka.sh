#!/bin/bash
gproject=$1
allcontigs=$2
refstrains=$3
outdir=$4
seqcentre=$5
if [ -z $seqcentre ] ; then 
  $seqcentre='XXX' 
fi

prokkadir=$(dirname `which prokka`)/../

echo "### assembly: $gproject; contig files from: ${allcontigs}/"
head -n1 $refstrains
taxo=(`grep $gproject $refstrains`)
echo ${taxo[@]}
genus=${taxo[1]} ; species=${taxo[2]%*.} ; strain=${taxo[3]} ; taxid=${taxo[4]} ; loctagprefix=${taxo[5]}
mkdir -p $outdir
if [ -e $prokkadir/db/genus/$genus ] ; then
 usegenus="--usegenus"
fi
prokkaopts="
--outdir $outdir --prefix ${genus}_${species}_${strain} --force 
--addgenes --locustag ${loctagprefix} --compliant --centre ${seqcentre} ${usegenus}
--genus ${genus} --species ${species} --strain '${strain}'
 --kingdom Bacteria --gcode 11 --cpu 4"
echo "#call: prokka $prokkaopts ${allcontigs}"
prokka $prokkaopts ${allcontigs}
date
