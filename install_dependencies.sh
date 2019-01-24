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
  if [ $? != 0 ]; then echo "ERROR: $1" ; exit 1 ; fi
}

if [ -z $1 ] ; then
  SOFTWARE=$PWD
else
  SOFTWARE=$(readlink -f $1)
fi
if [ -z $2 ] ; then
  BINS=${HOME}/bin
else
  BINS=$(readlink -f $2)
  PATH=${PATH}:${BINS}
fi


echo "Installation of Pantagruel and dependencies: ..."
echo ""

mkdir -p ${SOFTWARE} ${BINS}

cd ${SOFTWARE}

# Install Pantagruel pipeline and specific phylogenetic modules 
echo "get/update git repositories for Pantagruel pipeline and specific phylogenetic modules"
for repo in tree2 pantagruel ; do
  if [ -d ${SOFTWARE}/${repo} ] ; then
    cd ${repo} && git pull && cd -
  else
    git clone https://github.com/flass/${repo}
  fi
  if [ $? != 0 ] ; then
    echo "ERROR: Unable to create/update git repository ${SOFTWARE}/${repo}"
    exit 1
  fi
done
echo ""

# basic dependencies, libs and standalone software, R and packages, Python and packages
echo "get/update required Debian packages"
deppackages="git build-essential cmake gcc g++ linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes r-base-core r-recommended r-cran-ape r-cran-phytools r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite r-cran-igraph r-cran-getopt python python-scipy python-numpy python-biopython python-igraph cython mpi-default-bin mpi-default-dev mrbayes-mpi docker.io python-pip openjdk-8-jdk openjdk-8-jre"
sudo apt install $deppackages
if [ $? != 0 ] ; then
  echo "ERROR: could not install all required Debian packages:"
  apt list $deppackages
  exit 1
fi
echo ""

# install Python package BCBio.GFF
sudo pip install bcbio-gff

# install R packages not supported by Debian
if [[ -z "$(Rscript -e 'print(installed.packages()[,1:2])' | grep topGO | cut -d' ' -f1)" ]] ; then
sudo R --vanilla <<EOF
source("https://bioconductor.org/biocLite.R")
biocLite("topGO")
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
touch ${HOME}/.bash_profile
if [[ -z "$(grep PATH ${HOME}/.bash_profile | grep linuxbrew)" ]] ; then
  echo 'export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"' >> ${HOME}/.bash_profile
fi
if [[ -z "$(grep MANPATH ${HOME}/.bash_profile | grep linuxbrew)" ]] ; then
  echo 'export MANPATH="/home/linuxbrew/.linuxbrew/share/man:$MANPATH"' >> ${HOME}/.bash_profile
fi
if [[ -z "$(grep INFOPATH ${HOME}/.bash_profile | grep linuxbrew)" ]] ; then
  echo 'export INFOPATH="/home/linuxbrew/.linuxbrew/share/info:$INFOPATH"' >> ${HOME}/.bash_profile
fi
if [[ -z "$(export | grep PATH | grep linuxbrew)" ]] ; then
  export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"
fi

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
  ln -s ${SOFTWARE}/pal2nal.v14/pal2nal.pl ${BINS}/
  checkexec "Could not install pal2nal.pl"
fi
echo ""

# fetch MAD program
echo "get MAD"
if [ ! -x ${SOFTWARE}/mad/mad ] ; then
  wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
  unzip ${SOFTWARE}/mad2-2.zip
  chmod +x ${SOFTWARE}/mad/mad
  ln -s ${SOFTWARE}/mad/mad ${BINS}/
  checkexec "Could not install MAD"
fi
echo ""

# set up Docker group (with root-equivalent permissions) and add main user to it
# !!! makes the system less secure: OK within a dedicated virtual machine but to avoid on a server or desktop (or VM with other use)
if [[ -z "$(grep docker /etc/group)" ]] ; then
  sudo groupadd docker
fi
if [[ -z "$(grep docker /etc/group | grep ${USER})" ]] ; then
  sudo usermod -aG docker $USER
  newgrp docker
  echo ""
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
if [[ -z "$(which ALEml)" ]] ; then
  alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"
  alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"
  alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"
fi
if [[ -z "$(grep 'alias ALEml' ${HOME}/.bashrc)" ]] ; then
  echo 'alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"' >> ${HOME}/.bashrc
  echo 'alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"' >> ${HOME}/.bashrc
  echo 'alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"' >> ${HOME}/.bashrc
  checkexec "Could not store command aliases in ~/.bashrc"
fi

# add Python modules to PYTHONPATH
if [[ -z "$(grep PYTHONPATH ${HOME}/.bashrc | grep tree2 )" || -z "$(grep PYTHONPATH ${HOME}/.bashrc | grep pantagruel/python_libs )" ]] ; then
  echo 'export PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/tree2:${SOFTWARE}/pantagruel/python_libs' >> ${HOME}/.bashrc
  checkexec "Could not store PYTHONPATH in .bashrc"
fi
if [[ -z "$(export | grep PYTHONPATH | grep tree2 )" || -z "$(export | grep PYTHONPATH | grep pantagruel/python_libs )" ]] ; then
  export PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/tree2:${SOFTWARE}/pantagruel/python_libs
fi
echo ""

# install Interproscan
echo "get Interproscan"
ipversion="5.32-71.0"
if [[ ! -x ${SOFTWARE}/interproscan-${ipversion}/interproscan.sh ]] ; then
  ipsourceftprep="ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${ipversion}/"
  ipsourcefile="interproscan-${ipversion}-64-bit.tar.gz"
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
  ${SOFTWARE}/interproscan-${ipversion}/interproscan.sh -i ${SOFTWARE}/interproscan-${ipversion}/test_proteins.fasta -f tsv
  checkexec "Interproscan test was not successful"
else
  echo "found InterProScan executable:"
  ls ${SOFTWARE}/interproscan-${ipversion}/interproscan.sh
fi


if [[ -z "$(grep 'PATH=$PATH:$HOME/bin' ${HOME}/.bashrc)" ]] ; then
  echo 'export PATH=$PATH:$HOME/bin' >> ${HOME}/.bashrc
fi
if [[ -z "$(export | grep ' PATH=' | grep $HOME/bin)" ]] ; then
  export PATH=$PATH:$HOME/bin
fi

mkdir -p ${HOME}/bin/
rm -f ${HOME}/bin/pantagruel
ln -s ${SOFTWARE}/pantagruel/pantagruel ${HOME}/bin/

echo "Installation of Pantagruel and dependencies: complete"
