#!/bin/bash
gproject=$1
allcontigs=$2
refstrains=$3
outdir=$4
seqcentre=$5
if [ -z $seqcentre ] ; then 
  $seqcentre='XXX' 
fi

letters='abcdefghijklmnopqrstuvwxyz'

prokkadir=$(dirname `which prokka`)/../

echo "### assembly: $gproject; contig files from: ${allcontigs}/"
head -n1 $refstrains
taxo=(`grep $gproject $refstrains`)
echo ${taxo[@]}
genus=${taxo[1]} ; species=${taxo[2]%*.} ; strain=${taxo[3]} ; taxid=${taxo[4]} ; loctagprefix=${taxo[5]}
nrepli=0
nplas=0
for fastarepli in `ls -S $allcontigs/*.fasta` ; do
repli=`basename $fastarepli`
plasmopt=""
let "nrepli = $nrepli + 1"
replistatus=`echo $repli | cut -d'.' -f2`
if [[ $replistatus=='circularised' || $replistatus=='complete' ]] ; then
 let "nplas = $nplas + 1"
 if [ $nrepli -gt 0 ] ; then
  plet=${letters:$nplas:1}
  replname="p${loctag}-${plet}"
  plasmopt="--plasmid $replname"
  loctag="${loctagprefix}_p${plet}${nrepli}"
 else
  replname="chr."
  loctag="${loctagprefix}_chr${nrepli}"
 fi
else
 loctag="${loctagprefix}_${nrepli}"
fi
replitag=`echo ${repli%%.*} | sed -e "s/\(${strain}[\._]\)*\(.\+\)/\2/g"`
echo "## contig #$nrepli"
echo "${genus}_${species}_${strain}.${replitag} ${replname}"
date
mkdir -p $outdir
if [ -e $prokkadir/db/genus/$genus ] ; then
 usegenus="--usegenus"
fi
prokkaopts="
--outdir $outdir --prefix ${genus}_${species}_${strain}.${replitag} --force 
--addgenes --locustag ${loctag} --compliant --centre ${seqcentre} ${usegenus}
--genus ${genus} --species ${species} --strain '${strain}' ${plasmopt}
 --kingdom Bacteria --gcode 11 --cpu 4"
echo "#call: prokka $prokkaopts ${allcontigs}/${repli}"
prokka $prokkaopts ${allcontigs}/${repli}
date
done
