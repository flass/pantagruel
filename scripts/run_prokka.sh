#!/bin/bash
gproject=$1
allcontigs=$2
refstrains=$3
outdir=$4
seqcentre=$5
if [ -z $seqcentre ] ; then 
  $seqcentre='XXX' 
fi

echo "### assembly: $gproject; contig files from: ${allcontigs}/"
head -n1 $refstrains

taxo=($(grep -P "^${gproject}\t" ${refstrains} | sed -e 's/ /@#=/g')) # escapes the whitespaces
# restores the whitespaces or replace them with underscores when deemsd unsafe
genus="$(echo ${taxo[1]} | sed -e 's/@#=/ /g')"
species="$(echo ${taxo[2]%*.} | sed -e 's/@#=/ /g')"
strain="$(echo ${taxo[3]} | sed -e 's/@#=/ /g')"
safestrain="$(echo ${taxo[3]} | sed -e 's/@#=/_/g')"
if [ "${strain}" != "${safestrain}" ] ; then
  echo -e "Error: it is not safe to have a strain name with a whitespace in it: '${strain}'.\nPlease edit your strain information file '${refstrains}' and re-run this step. Will exit now" >&2
  exit 1
fi
taxid="$(echo ${taxo[4]} | sed -e 's/@#=/_/g')"
loctagprefix="$(echo ${taxo[5]} | sed -e 's/@#=/_/g')"
echo "'${genus}' '${species}' '${strain}' '${taxid}' '${loctagprefix}'"

if [[ -z "${genus}" || -z "${species}" || -z "${strain}" || -z "${taxid}" || -z "${loctagprefix}" ]] ; then
  echo "Error: missing value in genus='$genus', species='$species', strain='$strain', taxid='$taxid', loctagprefix='$loctagprefix'" >&2
  echo "Cannot run Prokka, exit now." >&2
  exit 1
fi
mkdir -p ${outdir}

prokkabin=$(which prokka)
if [ ! -z "$(env | grep PROKKA_DATA_DIR | cut -d'=' -f2)" ] ; then
  prokkavers=$(${prokkabin} --version 2>&1 >/dev/null | tr ' ' '-')
  prokkablastdb=${PROKKA_DATA_DIR}/${prokkavers}/db/genus
  echo "The environment variable \${PROKKA_DATA_DIR} is set; using its value as the location of Prokka reference BLAST databases:" >&2
  ls ${prokkablastdb} > /dev/null
  if [ ! -d ${prokkablastdb} ] ; then
    echo "Error: the directory '${prokkablastdb}' is missing" >&2
	exit 1
  fi
else
  prokkadir=$(dirname $(dirname $(readlink -f ${prokkabin} | awk '{ print $NF }')))
  prokkablastdb=${prokkadir}/db/genus
fi

if [ -e ${prokkablastdb}/${refgenus} ] ; then
 usegenus="--usegenus"
 echo "Will use '${prokkablastdb}/${refgenus}' protein database for BLAST-based annotation"
 if [[ "${refgenus}" != "${genus}" ]] ; then
   # refgenus value is 'Reference', or different genus than that of the focal genome (hereafter 'Genus')
   if [[ -e ${prokkablastdb}/${genus} ]] ; then
     # substitute the 'Reference' files to the 'Genus' files for the time of this annotation
     echo "Temporarily move the resident genus database '${prokkablastdb}/${genus}' to '${prokkablastdb}/${genus}_bak' to substitute it with '${prokkablastdb}/${refgenus}' and make links:"
     for gdbfile in $(ls ${prokkablastdb}/${genus} ${prokkablastdb}/${genus}.*) ; do
       mv -f ${gdbfile} ${gdbfile/${genus}/${genus}_bak}
	   refdbfile=${prokkablastdb}/$(basename ${gdbfile/${genus}/${refgenus}})
       ln -s ${refdbfile} ${gdbfile}
       ls -l ${gdbfile}
     done
   else
     # substitute the 'Reference' files to the 'Genus' files for the time of this annotation
     echo "temporarily link the reference database '${prokkablastdb}/${refgenus}' to '${prokkablastdb}/${genus}' for detection by Prokka:"
     for refdbfile in $(ls ${prokkablastdb}/${refgenus} ${prokkablastdb}/${refgenus}.*) ; do
       gdbfile=${prokkablastdb}/$(basename ${refdbfile/${refgenus}/${genus}})
	   ln -s ${refdbfile} ${gdbfile}
	   ls -l ${gdbfile}
     done
   fi
 fi
elif [ -e ${prokkablastdb}/${genus} ] ; then
 usegenus="--usegenus"
else
 echo "Warning: could not find the genus-specific database file '${prokkablastdb}/${genus}' or the automaticly generated placeholder '${prokkablastdb}/${refgenus}'; will NOT use the \`--genus\` option i.e. will annotate as using the generic Bacteria db" >&2
fi

if [ ! -z "${ptgthreads}" ] ; then
  paraopt="--cpus ${ptgthreads}"
else
  paraopt="--cpus $(nproc)"
fi
prokkaopts="
--outdir ${outdir} --prefix ${genus}_${species}_${safestrain} --force 
--addgenes --locustag ${loctagprefix} --compliant --centre ${seqcentre} ${usegenus}
--genus ${genus} --species ${species} --strain ${strain}
 --kingdom Bacteria --gcode 11 ${paraopt}"
echo "#call: ${prokkabin} ${prokkaopts} ${allcontigs}"
${prokkabin} ${prokkaopts} ${allcontigs}
prokkaexit=${?}
date

# revert to the original 'Genus' database
for dbfile in $(ls ${prokkablastdb}/${genus} ${prokkablastdb}/${genus}.*) ; do
  if [ -L ${dbfile} ] ; then
    rm ${dbfile}
    echo "Removed link: ${dbfile}"
  fi
  if [ -e ${dbfile/${genus}/${genus}_bak} ] ; then
    echo "Restore the resident genus database file '${dbfile}' from '${dbfile/${genus}/${genus}_bak}'"
    mv -f ${dbfile/${genus}/${genus}_bak} ${dbfile}
  fi
done

if [ ${prokkaexit} -gt 0 ] ; then 
  echo "An error occurred during prokka run; see logs in $(ls ${outdir}/*.log)" >&2
  echo "DEBUG: try running the last command caled by Prokka, but keeping the STDERR:" >&2
  failedcmd=$(grep "Could not run command: " ${outdir}/*.log | sed -e 's#\[.\+\] Could not run command: \(.\+\) 2> /dev/null#\1#')
  echo "# ${failedcmd}" >&2
  eval "${failedcmd}"
  exit ${prokkaexit}
fi