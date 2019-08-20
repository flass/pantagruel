#! /bin/bash

if [ -z ${2} ] ; then
  echo "Usage: ${0} assembly_list ouput_dir [nb_threads, default=min(\`nproc\`,4)] [specific_assembly_file_glob_pattern]"
  exit 1
fi

asslist="${1}"
outdir="${2}"

if [ ! -z "${3}" ] ; then
  para="--parallel=${3}"
elif [ $(nproc) -gt 1 ] ; then
  if [ $(nproc) -le 4 ] ; then
    para="--parallel=$(nproc)"
  else
    para="--parallel=4"
  fi
else
  para=""
fi
if [ ! -z "${4}" ] ; then
  fpat="${4}"
else
  fpat=""
fi

### fetch genome assembly data from the NCBI Assembly database
user='anonymous'
pswd=${myemail}
host='ftp.ncbi.nlm.nih.gov'
openparam=" -u $user,$pswd $host"
cd ${outdir}/

for ass in $(cat ${asslist}) ; do
  assloc="/genomes/all/${ass:0:3}/${ass:4:3}/${ass:7:3}/${ass:10:3}"
  # first fetch the complete name of target assembly folder
  assdir=$(lftp $openparam -e "ls ${assloc}/${ass}* ; quit" | awk '{print $NF}')
  if [ -z ${assdir} ] ; then
    echo "could not find folder of accession '${ass}' on NCBI FTP when looking for folder matching '${assloc}/${ass}*' ; exit now"
    exit 1
  fi
  fullass=$(basename ${assdir})
  if [ ! -e ${outdir}/${fullass}/md5checksums.txt ] ; then
    echo -e "fetch ${assdir}\n"
    # then download
	if [ -z "${fpat}" ] ; then
	  # the target folder
	  lftp ${openparam} -e "set ftp:use-feat false ; mirror -X *_assembly_structure ${para} ${assdir} ; quit"
	else
	  # the specific target files
	  mkdir -p ${outdir}/${fullass}/
	  lftp ${openparam} -e "set ftp:use-feat false ; mget -O ${fullass} ${assdir}/${fpat} ; quit"
	fi
	ls -l ${outdir}/${fullass}/
  else
    echo "assembly ${fullass} already present in ${outdir}/"
  fi
  cd ${outdir}/${fullass}/
  md5sum -c md5checksums.txt
  if [ $? -gt 0 ] ; then
#    if [ ! -z "$(md5sum --quiet -c md5checksums.txt 2> /dev/null | grep -v 'assembly_structure')" ] ; then
      echo "Error: files in ${ncbiass}/${fullass}/ seem corrupted (not only about missing *assembly_structure/ files) ; exit now"
      exit 1
#    else
#      echo "Warning: some files are corrupted or missing in the *assembly_structure/ ; this is not important for Pantagruel though."
#    fi
  fi
  cd ${outdir}/
  echo -e "${ass}: done\n"
done
