#! /bin/bash

asslist="$1"
outdir="$2"
if [ !-z "$3" ] ; then
  para="--parallel=$3"
elif [ $(nproc) -gt 1 ] ; then
  if [ $(nproc) -le 4 ] ; then
    para="--parallel=$(nproc)"
  else
    para="--parallel=4"
  fi
else
  para=""
fi

### fetch REFERENCE genome assembly data from the NCBI Assembly database
user='anonymous'
pswd=${myemail}
host='ftp.ncbi.nlm.nih.gov'
openparam=" -u $user,$pswd $host"
cd ${outdir}/

for ass in $(cat ${asslist}) ; do
  assloc="/genomes/all/${ass:0:3}/${ass:4:3}/${ass:7:3}/${ass:10:3}"
  # first fetch the complete name of target assembly folders
  assdir=$(lftp $openparam -e "ls ${assloc}/${ass}* ; quit" | awk '{print $NF}')
  if [ -z ${assdir} ] ; then
    echo "could not find folder of accession '${ass}' on NCBI FTP when looking for folder matching '${assloc}/${ass}*' ; exit now"
    exit 1
  fi
  fullass=$(basename ${assdir})
  if [ ! -e ${outdir}/${fullass}/md5checksums.txt ] ; then
    echo -e "fetch ${assdir}\n"
    # then download it
    lftp ${openparam} -e "mirror ${para} ${assdir} ; quit"
    ls -l ${outdir}/${fullass}/
  else
    echo "assembly ${fullass} already present in ${outdir}/"
  fi
  cd ${outdir}/${fullass}/
  if [ -d ${fullass}_assembly_structure ] ; then
    sudo chown -R $USER:$USER ${fullass}_assembly_structure
    chmod -R +rX ${fullass}_assembly_structure
  fi
  md5sum -c md5checksums.txt
  if [ $? -gt 0 ] ; then
    if [ $(md5sum --quiet -c md5checksums.txt 2> /dev/null | grep -v 'assembly_structure') ] ; then
      echo "Error: files in ${ncbiass}/${fullass}/ seem corrupted (not only about missing *assembly_structure/ files) ; exit now"
      exit 1
    else
      echo "Warning: some files are correupted or missing in the *assembly_structure/ ; this is not important for Pantagruel though."
    fi
  fi
  cd ${outdir}/
  echo -e "${ass}: done\n"
done
