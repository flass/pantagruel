#!/usr/bin/env bash

#########################################################
## PANTAGRUEL:                                         ##
##             a pipeline for                          ##
##             phylogenetic reconciliation             ##
##             of a bacterial pangenome                ##
#########################################################
# Copyright: Florent Lassalle (f.lassalle@imperial.ac.uk), 01 October 2018.

CRANmirror="https://cran.ma.imperial.ac.uk/" # should be edited to match nearest CRAN mirror

checkexec (){
  if [ $? != 0 ]; then
    echo "ERROR: $1" 1>&2
    exit 1
  else
    if [ ! -z "$2" ] ; then
      echo "$2"
    fi
  fi
}

usage (){
	echo "Usage: $0 software_dir [bin_dir]"
	echo "  software_dir: target installation folder for pantagruel and other software"
	echo "  bin_dir: folder where to link relevant executable files, including the main `pantagruel` excutable (defaults to ~/bin/); the path to this folder will be permatnently added to your path (via editing your ~/.bashrc)"
}

if [ -z "$1" ] ; then
  echo "Error: you must at least specify the target installation dir for pantagruel and other software:"
  usage
  exit 1
else
  if [[ "$1" == '-h' || "$1" == '--help' ]] ; then
    usage
    exit 0
  else
    SOFTWARE=$(readlink -f "$1")
  fi
fi
if [ -z "$2" ] ; then
  BINS=${HOME}/bin
else
  BINS=$(readlink -f "$2")
  PATH=${PATH}:${BINS}
fi


echo "Installation of Pantagruel and dependencies: ..."
echo ""

mkdir -p ${SOFTWARE}/ ${BINS}/

cd ${SOFTWARE}/

# Install Pantagruel pipeline and specific phylogenetic modules 
echo "get/update git repositories for Pantagruel pipeline"
if [ -d ${SOFTWARE}/pantagruel ] ; then
  cd ${SOFTWARE}/pantagruel && git pull && cd -
else
  git clone https://github.com/flass/pantagruel
fi
if [ $? != 0 ] ; then
  echo "ERROR: Unable to create/update git repository ${SOFTWARE}/pantagruel"
  exit 1
fi
cd ${SOFTWARE}/pantagruel && git submodule init && git submodule update && cd -
if [ $? != 0 ] ; then
  echo "ERROR: Unable to create/update git submodules"
  exit 1
fi
echo ""

#~ for repo in tree2 pantagruel ; do
  #~ if [ -d ${SOFTWARE}/${repo} ] ; then
    #~ cd ${repo} && git pull && cd -
  #~ else
    #~ git clone https://github.com/flass/${repo}
  #~ fi
  #~ if [ $? != 0 ] ; then
    #~ echo "ERROR: Unable to create/update git repository ${SOFTWARE}/${repo}"
    #~ exit 1
  #~ fi
#~ done
#~ echo ""

# basic dependencies, libs and standalone software, R and packages, Python and packages
echo "get/update required Debian packages"
deppackages="git build-essential cmake gcc g++ linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes r-base-core r-recommended r-cran-ape r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite r-cran-igraph r-cran-getopt sqlite3 sqlite3-doc libmagick++-dev python python-scipy python-numpy python-biopython python-biopython-sql python-igraph cython bioperl mpi-default-bin mpi-default-dev mrbayes-mpi docker.io python-pip openjdk-8-jdk openjdk-8-jre cd-hit"
sudo apt install $deppackages
if [ $? != 0 ] ; then
  echo "ERROR: could not install all required Debian packages:"
  apt list $deppackages
  exit 1
fi
echo ""

# install Python package BCBio.GFF
sudo -H pip install bcbio-gff

# install R packages not supported by Debian
if [[ -z "$(Rscript -e 'print(installed.packages()[,1:2])' | grep topGO | cut -d' ' -f1)" ]] ; then
sudo R --vanilla <<EOF
source("https://bioconductor.org/biocLite.R")
biocLite("topGO")
EOF
fi
if [[ -z "$(Rscript -e 'print(installed.packages()[,1:2])' | grep phytools | cut -d' ' -f1)" ]] ; then
sudo R --vanilla <<EOF
install.packages('phytools', repos='${CRANmirror}')
EOF
fi
if [[ -z "$(Rscript -e 'print(installed.packages()[,1:2])' | grep pvclust | cut -d' ' -f1)" ]] ; then
sudo R --vanilla <<EOF
install.packages('pvclust', repos='${CRANmirror}')
EOF
fi
echo "library('topGO')" | R --vanilla &> /dev/null
checkexec "Could not install R package 'topGO'"
echo "library('pvclust')" | R --vanilla &> /dev/null
checkexec "Could not install R package 'pvclust'"
echo ""

# Configure Linuxbrew
if [ ! -e ${HOME}/.bash_profile ] ; then
  # because when .bash_profile exists, it superseeds .profile, which won't be loaded automatically
  # thus if creating .bash_profile, restaure the loading of .profile
  echo '. "$HOME/.profile"' > ${HOME}/.bash_profile
fi
if [[ -z "$(grep PATH ${HOME}/.bash_profile | grep linuxbrew)" ]] ; then
  echo 'export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"' >> ${HOME}/.bash_profile
  editedprofile=true
fi
if [[ -z "$(grep MANPATH ${HOME}/.bash_profile | grep linuxbrew)" ]] ; then
  echo 'export MANPATH="/home/linuxbrew/.linuxbrew/share/man:$MANPATH"' >> ${HOME}/.bash_profile
  editedprofile=true
fi
if [[ -z "$(grep INFOPATH ${HOME}/.bash_profile | grep linuxbrew)" ]] ; then
  echo 'export INFOPATH="/home/linuxbrew/.linuxbrew/share/info:$INFOPATH"' >> ${HOME}/.bash_profile
  editedprofile=true
fi
source ${HOME}/.bash_profile

# install Prokka using brew
if [[ -z "$(brew search prokka | grep prokka)" ]] ; then
  brew update
  brew doctor
  brew install brewsci/bio/prokka
  checkexec "Could not install Prokka using Brew"
else
  echo "found Prokka already installed with Brew:"
  brew search prokka
fi
echo ""

# install MMSeqs using brew
if [[ -z "$(brew search mmseqs2 | grep mmseqs2)" ]] ; then
  brew doctor
  brew install mmseqs2
  checkexec "Could not install MMSeqs using Brew"
else
  echo "found mmseqs2 already installed with Brew:"
  brew search mmseqs2
fi
echo ""

# fetch Pal2Nal script
echo "get pal2nal"
if [ ! -x ${SOFTWARE}/pal2nal.v14/pal2nal.pl ] ; then
  wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
  tar -xzf ${SOFTWARE}/pal2nal.v14.tar.gz
  chmod +x ${SOFTWARE}/pal2nal.v14/pal2nal.pl
  checkexec "Could not install pal2nal.pl" "Succesfully installed pal2nal.pl"
fi
if [ ! -x ${BINS}/pal2nal.pl ] ; then
  ln -s ${SOFTWARE}/pal2nal.v14/pal2nal.pl ${BINS}/
  checkexec "Could not link pal2nal.pl executable to ${BINS}/" "Succesfully linked pal2nal.pl executable to ${BINS}/"
fi
echo ""

# fetch MAD program
echo "get MAD"
if [ ! -x ${SOFTWARE}/mad/mad ] ; then
  wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
  unzip ${SOFTWARE}/mad2-2.zip
  chmod +x ${SOFTWARE}/mad/mad
  checkexec "Could not install MAD" "Succesfully installed MAD"
fi
if [ ! -x ${BINS}/mad ] ; then
  ln -s ${SOFTWARE}/mad/mad ${BINS}/
  checkexec "Could not link MAD executable to ${BINS}/" "Succesfully linked MAD executable to ${BINS}/"
fi
echo ""

# set up Docker group (with root-equivalent permissions) and add main user to it
# !!! makes the system less secure: OK within a dedicated virtual machine but to avoid on a server or desktop (or VM with other use)
if [[ -z "$(grep docker /etc/group)" ]] ; then
  sudo groupadd docker
fi
if [[ -z "$(grep docker /etc/group | grep ${USER})" ]] ; then
  sudo usermod -aG docker $USER
  newgrp docker << EOF
EOF
fi
checkexec "Could not set group 'docker' or let user '$USER' join it"


# install ALE using Docker --- OK within a virtual machine
if [[ -z "$(sudo docker images | grep boussau/alesuite)" ]] ; then
  sudo docker pull boussau/alesuite
  checkexec "Could not install ALE suite using Docker"
else
  echo "found ALE suite already installed with Docker:"
  sudo docker images | grep boussau/alesuite
fi
if [[ -z "$(grep 'alias ALEml=' ${HOME}/.bashrc)" ]] ; then
  echo 'alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"' >> ${HOME}/.bashrc
  echo 'alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"' >> ${HOME}/.bashrc
  echo 'alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"' >> ${HOME}/.bashrc
  checkexec "Could not store command aliases in ~/.bashrc"
  editedrc=true
fi

# add Python modules to PYTHONPATH
if [[ -z "$(grep 'export PYTHONPATH=' ${HOME}/.bashrc | grep ${SOFTWARE}/pantagruel/python_libs )" ]] ; then
  echo "export PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/pantagruel/python_libs" >> ${HOME}/.bashrc
  checkexec "Could not store PYTHONPATH in .bashrc"
  editedrc=true
fi
echo ""

# install Interproscan
iphost="ftp://ftp.ebi.ac.uk"
iploc="pub/software/unix/iprscan/5/"

currIPversion=$(interproscan --version 2> /dev/null | head -n 1 | sed -e 's/InterProScan version //')
lastIPversion=$(lftp -c "open ${iphost} ; ls -tr ${iploc} ; quit" | tail -n 1 | awk '{print $NF}')
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
    #~ export PATH=$(dirname ${goodjava}):$PATH
    #~ ${SOFTWARE}/interproscan-${currIPversion}/interproscan.sh
   #~ checkexec "Interproscan test was not successful"
   #~ if [[ -z "$(grep JAVA4INTERPROSCAN ${HOME}/.bashrc )" ]] ; then
     #~ echo "export JAVA4INTERPROSCAN=${goodjava}" >> ${HOME}/.bashrc
     #~ checkexec "Could not store JAVA4INTERPROSCAN in .bashrc"
     #~ editedrc=true
   #~ fi
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
  echo 'export PATH=$PATH:${BINS}' >> ${HOME}/.bashrc
  editedrc=true
fi

rm -f ${BINS}/pantagruel
ln -s ${SOFTWARE}/pantagruel/pantagruel ${BINS}/
checkexec "Failed to link pantagruel executable to ${BINS}/" "Succesfully linked pantagruel executable to ${BINS}/"

echo "Installation of Pantagruel and dependencies: complete"

#~ if [[ -z "$(export | grep PATH | grep linuxbrew)" ]] ; then
  #~ export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"
#~ fi
#~ if [[ -z "$(which ALEml)" ]] ; then
  #~ alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"
  #~ alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"
  #~ alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"
#~ fi
#~ if [[ -z "$(export | grep PYTHONPATH | grep tree2 )" || -z "$(export | grep PYTHONPATH | grep pantagruel/python_libs )" ]] ; then
  #~ export PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/tree2:${SOFTWARE}/pantagruel/python_libs
#~ fi
#~ if [[ -z "$(export | grep ' PATH=' | grep ${BINS})" ]] ; then
  #~ export PATH=$PATH:${BINS}
#~ fi

if [ "${editedprofile}" == 'true' ] ; then
  echo "the file '~/.bash_profile' was edited to set environment variables; please reload it with \`source ~/.bash_profile\` or reset your session for changes to take effect" 
fi
if [ "${editedrc}" == 'true' ] ; then
  echo "the file '~/.bashrc' was edited to set environment variables; please reload it with \`source ~/.bashrc\` or reset your session for changes to take effect" 
fi
