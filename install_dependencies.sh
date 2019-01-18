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

mkdir -p ${SOFTWARE} ${BINS}

cd ${SOFTWARE}

# Install Pantagruel pipeline and specific phylogenetic modules 
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
# basic dependencies, libs and standalone software, R and packages, Python and packages
deppackages="git cmake gcc g++ linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes r-base-core r-recommended r-cran-ape r-cran-phytools r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite r-cran-igraph r-cran-getopt python python-scipy python-numpy python-biopython python-igraph cython mpi-default-bin mpi-default-dev mrbayes-mpi docker.io"
apt install $deppackages
if [ $? != 0 ] ; then
  echo "ERROR: could not install all required Debian packages:"
  apt list $deppackages
  exit 1
fi

# install R packages not supported by Debian
R --vanilla <<EOF
install.packages('pvclust', repos='${CRANmirror}')
source("https://bioconductor.org/biocLite.R")
biocLite("topGO")
EOF
echo "library('topGO')" | R --vanilla
checkexec "Could not install R package 'topGO'"
echo "library('pvclust')" | R --vanilla
checkexec "Could not install R package 'pvclust'"

# set up Docker group (with root-equivalent permissions) and add main user to it
# !!! makes the system less secure: OK within a dedicated virtual machine but to avoid on a server or desktop (or VM with other use)

groupadd docker
usermod -aG docker $USER
newgrp docker
checkexec "Could not set group 'docker' or let user '$USER' join it"

# install MMSeqs using brew
brew install https://raw.githubusercontent.com/soedinglab/mmseqs2/master/Formula/mmseqs2.rb --HEAD
checkexec "Could not install MMSeqs using Brew"

# fetch Pal2Nal script
wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -xzf pal2nal.v14.tar.gz
chmod +x ${SOFTWARE}/pal2nal.v14/pal2nal.pl
ln -s ${SOFTWARE}/pal2nal.v14/pal2nal.pl ${BINS}/
checkexec "Could not install pal2nal.pl"

# fetch MAD program
wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
mkdir -p ${SOFTWARE}/MAD/
tar -xzf ../mad2-2.zip -d ${SOFTWARE}/MAD/
chmod +x ${SOFTWARE}/MAD/mad
ln -s ${SOFTWARE}/MAD/mad
checkexec "Could not install MAD"

# install ALE using Docker --- OK within a virtual machine
docker pull boussau/alesuite
checkexec "Could not install ALE using Docker"
alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"
alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"
alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"
echo 'alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"' >> ${HOME}/.bashrc
echo 'alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"' >> ${HOME}/.bashrc
echo 'alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"' >> ${HOME}/.bashrc
checkexec "Could not store command aliases in .bashrc"

# add Python modules to PYTHONPATH
PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/tree2:${SOFTWARE}/pantagruel/python_libs
echo 'PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/tree2:${SOFTWARE}/pantagruel/python_libs' >> ${HOME}/.bashrc
checkexec "Could not store PYTHONPATH in .bashrc"

# install Interproscan
ipversion="5.32-71.0"
ipsourceftprep="ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${ipversion}/"
ipsourcefile="interproscan-${ipversion}-64-bit.tar.gz"
wget ${ipsourceftprep}/${ipsourcefile}
wget ${ipsourceftprep}/${ipsourcefile}.md5
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
tar -pxvzf ${ipsourcefile}
checkexec "Could not uncompress Interproscan successfully"
interproscan-${ipversion}/interproscan.sh -i test_proteins.fasta -f tsv
checkexec "Interproscan test was not successful"

# get the GO term db - not necessary so far
# wget http://archive.geneontology.org/latest-termdb/go_daily-termdb-data.gz
