#!/usr/bin/bash
SOFTWARE=$1
if [ -z ${SOFTWARE} ] ; then
  SOFTWARE=$PWD
fi
BINS=$2
if [ -z ${BINS} ] ; then
  BINS=${HOME}/bin
else
  PATH=${PATH}:${BINS}
fi

mkdir -p ${BINS}

cd ${SOFTWARE}


# basic dependencies, libs and standalone software
sudo apt install git cmake gcc g++ linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes
# R and packages
sudo apt install r-base-core r-recommended r-cran-ape r-cran-phytools r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite r-cran-igraph r-cran-getopt
# Python and packages
sudo apt install python python-scipy python-numpy python-biopython python-igraph cython

# if using MPI for MrBayes (recommended)
sudo apt install mpi-default-bin mpi-default-dev mrbayes-mpi

# if using Docker to install ALE
# --- OK within a virtual machine
# !!! to avoid on a server or desktop due to the need to grant root-equivalent right to main user
sudo apt install docker.io
sudo groupadd docker
sudo usermod -aG docker $USER
sudo newgrp docker

# install MMSeqs using brew
brew install https://raw.githubusercontent.com/soedinglab/mmseqs2/master/Formula/mmseqs2.rb --HEAD

# fetch Pal2Nal script
wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -xzf pal2nal.v14.tar.gz
chmod +x ${SOFTWARE}/pal2nal.v14/pal2nal.pl
ln -s ${SOFTWARE}/pal2nal.v14/pal2nal.pl ${BINS}/

# fetch MAD program
wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
mkdir -p ${SOFTWARE}/MAD/
tar -xzf ../mad2-2.zip -d ${SOFTWARE}/MAD/
chmod +x ${SOFTWARE}/MAD/mad
ln -s ${SOFTWARE}/MAD/mad

# install ALE using Docker --- OK within a virtual machine
docker pull boussau/alesuite
alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"
alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"
alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"
echo 'alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"' >> ${HOME}/.bashrc
echo 'alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"' >> ${HOME}/.bashrc
echo 'alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"' >> ${HOME}/.bashrc

# Install Pantagruel pipeline and specific phylogenetic modules 
cd ${SOFTWARE}/
git clone https://github.com/flass/tree2
git clone https://github.com/flass/pantagruel
```
PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/tree2:${SOFTWARE}/pantagruel/python_libs
echo 'PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/tree2:${SOFTWARE}/pantagruel/python_libs' >> ${HOME}/.bashrc

