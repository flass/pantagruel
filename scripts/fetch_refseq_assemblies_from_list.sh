#! /bin/bash

asslist="$1"
outdir="$2"

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
    echo "# lftp ${openparam} -e \"mirror ${assdir} ; quit\""
    lftp ${openparam} -e "mirror ${assdir} ; quit"
    ls -l ${outdir}/${fullass}/
  else
    echo "assembly ${fullass} already present in ${ncbiass}/"
  fi
  cd ${outdir}/${fullass}/
  if [ -d ${fullass}_assembly_structure ] ; then
    sudo chown -R $USER:$USER ${fullass}_assembly_structure
    chmod -R +rX ${fullass}_assembly_structure
  fi
  md5sum -c md5checksums.txt
  if [ $? -gt 0 ] ; then
    echo "Error: files in ${ncbiass}/${fullass}/ seem corrupted"
  fi
  cd ${ncbiass}/
  echo -e "${ass}: done\n"
done
