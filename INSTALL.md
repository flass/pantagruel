## Installing Pantagruel and its dependencies

Under a Debian environment (e.g. Ubuntu), please follow the indications below.  
In the case of setting-up a virtual machine, the script [install_dependencies.sh](https://github.com/flass/pantagruel/blob/master/install_dependencies.sh) can do that work for you.  

Assuming you want to create a folder named `pantagruel_pipeline/` in the current working directory and install the whole software pipeline (the *pantagruel* package and its dependencies) in it, you should first do:
```sh
mkdir ./pantagruel_pipeline/
```
Then, you have to get the pantagruel scripts by dowloading the [archive of the last version on Github](https://github.com/flass/pantagruel/archive/master.zip) or use `git` to synchronize the repository (recomended for easier software update, especially during *Pantagruel* development phase!).
```sh
cd pantagruel_pipeline/
git clone https://github.com/flass/pantagruel.git
```
Finally, you can run the installation script:  
!!! Note that for this step, you (i.e. your linux user) need to have the sudo rights (be part of the sudo group, check with `grep sudo /etc/group`); however, DO NOT execute the installtion script with `sudo` (which would make all software owned by the root user).  
```sh
cd .. # retrun to parent folder
pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline/
```  

Optionally, you can also specify the folder where relevant executable files, including the main `pantagruel` excutable, will be linked (it defaults to `~/bin/`, i.e. is user-specific). The path to this folder will be permatnently added to your path (via editing your ~/.bashrc).  
```sh
pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline/ pantagruel_pipeline/bin/
```
________

Otherwise, you can manually install the software following the indications below:

### basic dependencies, libraries and standalone software
```sh
sudo apt install git cmake gcc g++ libmagick++-dev linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes openjdk-8-jdk openjdk-8-jre cd-hit
```
### R and packages
```sh
sudo apt install r-base-core r-recommended r-cran-ape r-cran-phytools r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite r-cran-igraph r-cran-getopt
```
NB: debian package `r-cran-phytools` is only avalable from Ubunti 18.04; if running an earlier version, please
For optional packages `topGO` and `pvclust`, you can open the R interpreter using `R --vanilla` and typing:
```R
CRANmirror='https://cran.ma.imperial.ac.uk/'
source("https://bioconductor.org/biocLite.R")
biocLite("topGO")
install.packages('pvclust', repos=CRANmirror)
# and in case of not being able to install phytools via apt:
install.packages('phytools', repos=CRANmirror)
```

### Python and packages
```sh
sudo apt install python python-scipy python-numpy python-biopython python-igraph cython python-pip
```

### if using MPI for MrBayes (recommended)
```sh
sudo apt install mpi-default-bin mpi-default-dev mrbayes-mpi
```

### Fetch Pantagruel pipeline and specific phylogenetic modules
```sh
cd ${SOFTWARE}/
git clone https://github.com/flass/tree2
git clone https://github.com/flass/pantagruel
```
This commands need to be added to your .bashrc file for it to last beyond the current session:
```sh
PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/tree2:${SOFTWARE}/pantagruel/python_libs
```

### Install Prokka using brew
```sh
brew install brewsci/bio/prokka
```

### Install MMSeqs using brew
```sh
brew install mmseqs2
```

### Fetch Pal2Nal script
```sh
wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -xzf pal2nal.v14.tar.gz
chmod +x ${SOFTWARE}/pal2nal.v14/pal2nal.pl
ln -s ${SOFTWARE}/pal2nal.v14/pal2nal.pl ${BINS}/
```

### Fetch MAD program
```sh
wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
mkdir -p ${SOFTWARE}/MAD/
tar -xzf ../mad2-2.zip -d ${SOFTWARE}/MAD/
chmod +x ${SOFTWARE}/MAD/mad
ln -s ${SOFTWARE}/MAD/mad
```

### Using Docker to install ALE
This is recommended for use within a virtual machine with no other function.
It is however to avoid on a server or desktop due to the need to grant root-equivalent right to main user which leads to significant security breach exposure.
First install Docker:
```sh
sudo apt install docker.io
sudo groupadd docker
sudo usermod -aG docker $USER
sudo newgrp docker
```
Then install ALE:  
```sh
docker pull boussau/alesuite
```
Then set the following command aliases (add these lines to your .bashrc file for them to last beyond the current session):
```sh
alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"
alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"
alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"
```

### Installing ALE from source (s)
This is more cumbersome, but safer when installing this software on a computer that has other functions, especially on  al large production server. See [ALE's own installation page](https://github.com/ssolo/ALE/blob/master/INSTALL.md) for more details and updates.
```sh
sudo apt install libboost-dev libboost-serialization-dev libboost-mpi-dev
```
#### Bio++ libraries version >2.2.0 are required
The version 2.4.1 can be found as a Debian package on Ubuntu 18.4 LTS (Bionic Beaver). Previous Ubuntu versions, such as 16.4 LTS (Xenial Xerius) have version 2.1.0, in which case Bio++ have to be compiled from source (much heavier).  
Using Debian packages:  
```sh
sudo apt install libbpp-core-dev libbpp-phyl-dev libbpp-seq-dev libbpp-seq-omics-dev
```
OR from source:  
```sh
mkdir ${SOFTWARE}/bpp
cd ${SOFTWARE}/bpp
git clone https://github.com/BioPP/bpp-core
git clone https://github.com/BioPP/bpp-seq
git clone https://github.com/BioPP/bpp-phyl
mkdir bpp-core-build
mkdir bpp-phyl-build
mkdir bpp-seq-build
cd bpp-core-build/
cmake ../bpp-core -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
make
sudo make install
cd ../bpp-seq-build/
cmake ../bpp-seq -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
make
sudo make install
cd ../bpp-phyl-build/
cmake ../bpp-phyl -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
make
sudo make install
cd ../..
```
Finally install ALE:  
```sh
git clone https://github.com/ssolo/ALE
mkdir build
cd build
cmake ..
make
```
