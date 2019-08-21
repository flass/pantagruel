#! /bin/bash

if [ -z ${2} ] ; then
  echo "Usage: ${0} assembly_list ouput_dir [nb_threads, default=min(\`nproc\`,4)] [specific_assembly_file_suffix]"
  echo "   'specific_assembly_file_suffix' will be appended to the expected assembly file prefix, i.e. ASSEMBACESSION_ASSEMBNAME."
  echo "   The suffix can be a globable pattern."
  echo "   For instance, for the assembly with prefix 'GCA_000006745.1_ASM674v1':"
  echo "    - specifying the suffix '_genomic.fna.gz' will fetch only the file 'GCA_000006745.1_ASM674v1_genomic.fna.gz'."
  echo "    - specifying the suffix '*_genomic.fna.gz' will fetch several files: 'GCA_000006745.1_ASM674v1_genomic.fna.gz',"
  echo "      but also 'GCA_000006745.1_ASM674v1_cds_from_genomic.fna.gz' and 'GCA_000006745.1_ASM674v1_rna_from_genomic.fna.gz'."
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
  suf="${4}"
  # make a grep-freindly regexp out of pottential glob pattern
  sufpat="$(echo $suf | sed -e 's/*/.*/g')"
  suftxt="$(echo $suf | sed -e 's/*/_all/g')"
else
  suf=""
  sufpat=""
fi

### fetch genome assembly data from the NCBI Assembly database
user='anonymous'
pswd=${myemail}
host='ftp.ncbi.nlm.nih.gov'
openparam=" -u $user,$pswd $host"
cd ${outdir}/

fetchass (){
  echo -e "fetch ${assdir}\n"
  # then download
  if [ -z "${suf}" ] ; then
	# the target folder
    lftp ${openparam} -e "set ftp:use-feat false ; mirror -X *_assembly_structure ${para} ${assdir} ; quit"
  else
	# the specific target files
	mkdir -p ${outdir}/${fullass}/
	lftp ${openparam} -e "set ftp:use-feat false ; mget -O ${fullass} ${assdir}/${fullass}${suf} ; get -O ${fullass} ${assdir}/md5checksums.txt ; quit"
	grep ${fullass}${sufpat} ${outdir}/${fullass}/md5checksums.txt > ${outdir}/${fullass}/md5checksums${suftxt}.txt
  fi
  ls -l ${outdir}/${fullass}/
  cd ${outdir}/${fullass}/
  md5sum -c md5checksums${suftxt}.txt
  if [ $? -gt 0 ] ; then
    if [ ! -z "$(md5sum --quiet -c md5checksums${suftxt}.txt 2> /dev/null | grep -v 'assembly_structure')" ] ; then
      echo "Error: files in ${ncbiass}/${fullass}/ seem corrupted (not only about missing *assembly_structure/ files) ; exit now"
      exit 1
    else
      echo "Warning: some files are corrupted or missing in the *assembly_structure/ ; this is not important for Pantagruel though."
    fi
  fi
  cd ${outdir}/
}



for ass in $(cat ${asslist}) ; do
  assloc="/genomes/all/${ass:0:3}/${ass:4:3}/${ass:7:3}/${ass:10:3}"
  # first fetch the complete name of target assembly folder
  assdir=$(lftp $openparam -e "ls ${assloc}/${ass}* ; quit" | awk '{print $NF}')
  if [ -z ${assdir} ] ; then
    echo "could not find folder of accession '${ass}' on NCBI FTP when looking for folder matching '${assloc}/${ass}*' ; exit now"
    exit 1
  fi
  fullass=$(basename ${assdir})
  if [ -e ${outdir}/${fullass}/md5checksums.txt ] ; then
    cd ${outdir}/${fullass}/
    md5sum -c md5checksums.txt
	md5stat=${?}
    cd ${outdir}/
    if [ ${md5stat} -gt 0 ] ; then
	  echo "the folder of accession '${ass}' already exists locally, but the file set is not complete; try and repeat the transfer"
	  fetchass ${assdir}
    else
      echo "assembly ${fullass} already present in ${outdir}/"
    fi
  else
    fetchass ${assdir}
  fi
  echo -e "${ass}: done\n"
done
