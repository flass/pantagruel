#!/bin/bash
### script to install Interproscan 5 latest version

iphost="ftp://ftp.ebi.ac.uk"
iploc="pub/software/unix/iprscan/5/"

usage (){
	echo "Usage: $0 software_dir [bin_dir]"
	echo "  software_dir: target installation folder for pantagruel and other software"
	echo "  bin_dir: folder where to link the `interproscan` executable file (defaults to ~/bin/); the path to this folder will be permatnently added to your path (via editing your ~/.bashrc)"
}

if [ -z "${1}" ] ; then
  echo "Error: you must at least specify the target installation dir for pantagruel and other software:"
  usage
  exit 1
else
  if [[ "${1}" == '-h' || "${1}" == '--help' ]] ; then
    usage
    exit 0
  else
    SOFTWARE=$(readlink -f "${1}")
  fi
fi
shift

if [ ! -z "${2}" ] ; then
  BINS=$(readlink -f "${2}")
else
  BINS=${HOME}/bin
fi

PATH=${PATH}:${BINS}

echo "Installation of Interproscan"
echo ""

mkdir -p ${SOFTWARE}/ ${BINS}/

currIPversion=$(interproscan --version 2> /dev/null | head -n 1 | sed -e 's/InterProScan version //')
lastIPversion=$(lftp -c "open ${iphost} ; ls -tr ${iploc}/[0-9].[0-9][0-9]*-[0-9][0-9]*.[0-9] ; quit" | tail -n 1 | awk '{print $NF}')
if [[ -z "${currIPversion}" || "${currIPversion}" != "${lastIPversion}" ]] ; then
 echo "get Interproscan ${lastIPversion}"
 if [[ ! -x ${SOFTWARE}/interproscan-${lastIPversion}/interproscan.sh ]] ; then
   ipsourceftprep="${iphost}/${iploc}/${lastIPversion}/"
   ipsourcefile="interproscan-${lastIPversion}-64-bit.tar.gz"
   if [[ ! -e ${SOFTWARE}/${ipsourcefile} ]] ; then
     wget ${ipsourceftprep}/${ipsourcefile}
   fi
   if [[ ! -e ${SOFTWARE}/${ipsourcefile}.md5 ]] ; then
     wget ${ipsourceftprep}/${ipsourcefile}.md5
   fi
   dlok=$(md5sum -c ${ipsourcefile}.md5 | grep -o 'OK')
   ip=3
   # allow 3 attempts to downlaod the package, given its large size
   while [[ "$dlok" != 'OK' && $ip -gt 0 ]] ; do 
     ip=$(( $ip - 1 ))
     echo "dowload of ${ipsourcefile} failed; will retry ($ip times left)"
     wget ${ipsourceftprep}/${ipsourcefile}
     wget ${ipsourceftprep}/${ipsourcefile}.md5
     dlok=$(md5sum -c ${ipsourcefile}.md5 | grep -o 'OK')done
   done
   if [[ "$dlok" != 'OK' ]] ; then
     echo "ERROR: Could not dowload ${ipsourcefile}"
     exit 1
   fi
   tar -pxvzf ${SOFTWARE}/${ipsourcefile}
   checkexec "Could not uncompress Interproscan successfully"
 else
   echo "found InterProScan executable:"
   ls ${SOFTWARE}/interproscan-${lastIPversion}/interproscan.sh
 fi
 currIPversion=${lastIPversion}
 # make sure InterProscan uses the correct version of Java (1.8)
 wrongjava=$(${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh | grep -o 'Java version .* is required to run InterProScan.')
 if [ ! -z "${wrongjava}" ] ; then
   needjava=$(echo ${wrongjava} | sed -e 's/Java version \(.*\) is required to run InterProScan./\1/')
   currjava=$(readlink -f `which java`)
   goodjava=$(echo ${currjava} | sed -e "s/java-[0-9]\+-/java-${needjava}-/")
   echo "default version of Java is wrong: ${currjava}"
   echo "will find executable for desired version: ${needjava}"
   if [ ! -e ${goodjava} ] ; then
     goodjava=$(echo ${currjava} | sed -e "s/java-[0-9]\+-/java-${needjava}.0-/")
   fi
   if [ ! -e ${goodjava} ] ; then
     goodjava=$(ls -d $(dirname $(dirname $(dirname ${currjava})))/*${needjava}*/bin/java)
   fi
   if [ -z ${goodjava} ] ; then
     echo "ERROR: Could not find the required version of java for InterProScan:"
     ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh
     exit 1
   else
     echo "found version of java suitable for InterProScan: ${goodjava}"
   fi
   echo "will edit definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh"
   if [ ! -z "$(grep 'JAVA=\$(type -p java)' ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh)" ] ; then
     mv ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh.original && \
      sed -e "s#^JAVA=\$(type -p java)#JAVA=${goodjava}#g" ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh.original > ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh && \
      chmod +x ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh && \
      echo "fixed definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh (original script saved as '${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh.original')"
     checkexec "Could not edit definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh"
   else
     # the file was already modified, be more agressive with sed substitution
     mv ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh0 && \
      sed -e "s#^JAVA=.\+#JAVA=${goodjava}#g" ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh0 > ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh && \
      chmod +x ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh && rm ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh0 && \
      echo "fixed definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh"
     checkexec "Could not edit definition of JAVA in ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh"
   fi     
 fi
 ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh -i ${SOFTWARE}/interproscan-${currIPversion}/test_proteins.fasta -f tsv
 checkexec "Interproscan test was not successful"
 # link the exec file
 rm -f ${BINS}/interproscan
 ln -s ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh ${BINS}/interproscan
 checkexec "Failed to link Interproscan executable to ${BINS}/" "Succesfully linked Interproscan executable to ${BINS}/"
else
 echo "found up-to-date version of Interproscan at $(ls -l `which interproscan` | awk '{print $NF}')"
fi

if [[ -z "$(grep 'export PATH=' ${HOME}/.bashrc | grep ${BINS})" ]] ; then
  echo "export PATH=\$PATH:${BINS}" >> ${HOME}/.bashrc
  editedrc=true
fi

echo "Installation of InterProScan: complete"

if [ "${editedrc}" == 'true' ] ; then
  echo "the file '~/.bashrc' was edited to set environment variables; please reload it with \`source ~/.bashrc\` or reset your session for changes to take effect" 
fi
