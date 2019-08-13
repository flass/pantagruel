#!/bin/bash
gproject=$1
allcontigs=$2
refstrains=$3
outdir=$4
seqcentre=$5
if [ -z $seqcentre ] ; then 
  $seqcentre='XXX' 
fi

prokkablastdb=$(dirname $(dirname $(readlink -f `which prokka` | awk '{ print $NF }')))/db/genus

echo "### assembly: $gproject; contig files from: ${allcontigs}/"
head -n1 $refstrains
taxo=(`grep -P "^${gproject}\t"  ${refstrains}`)
echo ${taxo[@]}
genus=${taxo[1]} ; species=${taxo[2]%*.} ; strain=${taxo[3]} ; taxid=${taxo[4]} ; loctagprefix=${taxo[5]}
if [[ -z $genus || -z $species || -z $strain || -z $taxid || -z $loctagprefix ]] ; then
  echo "missing value in genus='$genus', species='$species', strain='$strain', taxid='$taxid', loctagprefix='$loctagprefix'"
  echo "cannot run Prokka, exit now."
  exit 1
fi
mkdir -p ${outdir}
if [ -e ${prokkablastdb}/${refgenus} ] ; then
 usegenus="--usegenus"
 if [ -e ${prokkablastdb}/${genus} ] ; then
   # substitute the 'Reference' files to the 'Genus' files for the time of this annotation
   echo "temporarily move the resident genus database '${prokkablastdb}/${genus}' to '${prokkablastdb}/${genus}_bak' to substitute '${prokkablastdb}/${refgenus}' to it"
   for gdbfile in $(ls ${prokkablastdb}/${genus} ${prokkablastdb}/${genus}.*) ; do
     mv -f ${gdbfile} ${gdbfile/${genus}/${genus}_bak}
     ln -s ${prokkablastdb}/$(basename ${gdbfile/${genus}/${refgenus}}) ${gdbfile}
   done
 else
   # substitute the 'Reference' files to the 'Genus' files for the time of this annotation
   for refdbfile in $(ls ${prokkablastdb}/${refgenus} ${prokkablastdb}/${refgenus}.*) ; do
     ln -s ${refdbfile} ${prokkablastdb}/$(basename ${refdbfile/${refgenus}/${genus}})
   done
 fi
 echo "made links:"
 ls -l ${gdbfile} ${gdbfile}.*
elif [ -e ${prokkablastdb}/${genus} ] ; then
 usegenus="--usegenus"
fi

prokkaopts="
--outdir $outdir --prefix ${genus}_${species}_${strain} --force 
--addgenes --locustag ${loctagprefix} --compliant --centre ${seqcentre} ${usegenus}
--genus ${genus} --species ${species} --strain '${strain}'
 --kingdom Bacteria --gcode 11"
echo "#call: prokka $prokkaopts ${allcontigs}"
prokka $prokkaopts ${allcontigs}
date

# revert to the original 'Genus' database
for dbfile in $(ls ${prokkablastdb}/${genus} ${prokkablastdb}/${genus}.*) ; do
  if [ -L ${dbfile} ] ; then
    rm ${dbfile}
    echo "removed link: ${dbfile}"
  fi
  if [ -e ${dbfile/${genus}/${genus}_bak} ] ; then
    echo "restore the resident genus database file '${dbfile}' from '${dbfile/${genus}/${genus}_bak}'"
    mv -f ${dbfile/${genus}/${genus}_bak} ${dbfile}
  fi
done
