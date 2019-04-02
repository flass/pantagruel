# Installing Pantagruel and its dependencies

## The automatic, easy way

Under a Debian environment (e.g. Ubuntu), you can automatically install all dependencies using the script [install_dependencies.sh](https://github.com/flass/pantagruel/blob/master/install_dependencies.sh), following the indications below:  

Assuming you want to create a folder named `pantagruel_pipeline/` in the current working directory and install the whole software pipeline (the *pantagruel* package and its dependencies) in it, you should first do:
```sh
mkdir ./pantagruel_pipeline/
```
Then, you have to get the pantagruel scripts by dowloading the [archive of the last version on Github](https://github.com/flass/pantagruel/archive/master.zip) or use `git` to synchronize the repository (recomended for easier software update, especially during *Pantagruel* development phase!).
```sh
cd pantagruel_pipeline/
git clone https://github.com/flass/pantagruel.git
```
Finally, you may run the installation script:  
!!! Note that for this step, you (i.e. your linux user) need to have the sudo rights (be part of the sudo group, check with `grep sudo /etc/group`); however, DO NOT execute the installtion script with `sudo` (which would make all software owned by the root user).  
```sh
cd .. # return to parent folder
pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline/
```  

Optionally, you can also specify the folder where relevant executable files, including the main `pantagruel` excutable, will be linked (it defaults to `~/bin/`, i.e. is user-specific). The path to this folder will be permatnently added to your path (via editing your ~/.bashrc).  
```sh
pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline/ pantagruel_pipeline/bin/
```
________

## The details if you are picky

If you are using a Debian environment but don't want to run the install script all at once (you may have your own reasons), here is a summary of what it does, so you can manually install each piece of software at your own discretion.  
These indications would also stand for installation on another type of Linux/Unix OS (e.g. Redhat, MacOS), but you would have to adapt the first part (using `apt` to install system packages on a Debian OS) to the relevant system (e.g. for a Red Hat OS, using `yum` instead of `apt`, and finding the matching packages).  
A special note clarifying deferent scenarios to install ALE can be found [here](https://github.com/flass/pantagruel/blob/master/doc/installing_ALE.md).

### Platform-dependent software installation: using Debian system packages

Most required software can simply be installed using the Debian system package manager `apt` (or any frontend like `synaptic`):    

#### basic dependencies, libraries and standalone software
```sh
sudo apt install git cmake gcc g++ libmagick++-dev sqlite3 sqlite3-doc linuxbrew-wrapper bioperl lftp clustalo raxml libhmsbeagle1v5 mrbayes openjdk-8-jdk openjdk-8-jre cd-hit
```
#### R and packages
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

#### Python and packages
```sh
# install Python 2.7, the Cython extension, the pip package manager and core packages
sudo apt install python python-scipy python-numpy python-biopython python-biopython-sql python-igraph cython python-pip
```

#### if using MPI for MrBayes (recommended)
```sh
sudo apt install mpi-default-bin mpi-default-dev mrbayes-mpi
```

Regarding the installation of ALE, you may need to install either the Docker daemon, or the Boost and Bio++ C++ libraries (see [installing_ALE](https://github.com/flass/pantagruel/blob/master/doc/installing_ALE.md)).

#### Install Docker:
```sh
sudo apt install docker.io
sudo groupadd docker
sudo usermod -aG docker $USER
sudo newgrp docker
```
#### Install Boost libraries
This is more cumbersome, but safer when installing this software on a computer that has other functions, especially on  al large production server. See [ALE's own installation page](https://github.com/ssolo/ALE/blob/master/INSTALL.md) for more details and updates.
```sh
sudo apt install libboost-dev libboost-serialization-dev libboost-mpi-dev
```
#### Install Bio++ libraries 
To compile ALE from source, Bio++ version >2.2.0 is required. The version 2.4.1 can be found as a Debian package on Ubuntu 18.4 LTS (Bionic Beaver). Previous Ubuntu versions, such as 16.4 LTS (Xenial Xerius) have version 2.1.0, in which case Bio++ have to be compiled from source (much heavier).  
```sh
sudo apt install libbpp-core-dev libbpp-phyl-dev libbpp-seq-dev libbpp-seq-omics-dev
```

### Platform-independent software installation

#### Fetch Pantagruel pipeline and specific phylogenetic modules
```sh
# assuming you want to install this in a folder named 'software':
cd software/
git clone https://github.com/flass/pantagruel
cd pantagruel/ && git submodule init && git submodule update && cd -
```
You then need to add the Python modules to your `PYTHONPATH`.
This commands need to be added to your `~/.bashrc` file for the module to be permanently available:
```sh
PYTHONPATH=${PYTHONPATH}:/full/path/to/software/pantagruel/python_libs
```

#### Install Prokka using brew
```sh
brew install brewsci/bio/prokka
```

#### Install MMSeqs using brew
```sh
brew install mmseqs2
```

#### Fetch Pal2Nal script
```sh
wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -xzf pal2nal.v14.tar.gz
chmod +x ${SOFTWARE}/pal2nal.v14/pal2nal.pl
ln -s ${SOFTWARE}/pal2nal.v14/pal2nal.pl ${BINS}/
```

#### Fetch MAD program
```sh
wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
mkdir -p ${SOFTWARE}/MAD/
tar -xzf ../mad2-2.zip -d ${SOFTWARE}/MAD/
chmod +x ${SOFTWARE}/MAD/mad
ln -s ${SOFTWARE}/MAD/mad
```

#### Install custom packages from PyPI
```sh
# install BCBio.GFF
sudo -H pip install bcbio-gff

# install bioscripts.convert
sudo -H pip install bioscripts.convert
```

#### Using Docker to install ALE
This is recommended for use within a virtual machine with no other function.
It is however to avoid on a server or desktop due to the need to grant root-equivalent right to main user which leads to significant security breach exposure.

```sh
docker pull boussau/alesuite
```
Then set the following command aliases (add these lines to your .bashrc file for them to last beyond the current session):
```sh
alias ALEobserve="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve"
alias ALEml="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml"
alias ALEml_undated="docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated"
```
______________________________________________

## List of software dependencies

If you feel like installing the depedencies your own way without following the above instructions, here is a list of the software on which Pantagruel depends:

### Required bioinformatic software
- **Prokka** for genome annotation
  (Install from  [source code](https://github.com/tseemann/prokka) or as brew package *prokka*)

- **MMseqs2/Linclust** for homologous sequence clustering  
  (Install from [source code](https://github.com/soedinglab/MMseqs2) or as brew package *mmseqs*)

- **Clustal Omega** for homologous sequence alignment  
  (Install from [source code](http://www.clustal.org/omega/) or as Debian package *clustalo*; version used and recommended: 1.2.1)  
  - \[ future development: consider using [FAMSA](http://sun.aei.polsl.pl/REFRESH/famsa) \]

- **PAL2NAL** for reverse tanslation of protein sequence alignments into CDS alignments  
  ([Perl source code](http://www.bork.embl.de/pal2nal/))

- **RAxML** for species tree and initial (full) gene tree estimation  
  (Install from [source code](https://github.com/stamatak/standard-RAxML) or as Debian package *raxml*; version used and recommended: 8.2.9)  
  - \[ future development: consider using RAxML-NG (Install from [source code](https://github.com/amkozlov/raxml-ng)) \]

- **MrBayes** for secondary estimation of (collapsed) gene trees  
  (Install from [source code](http://mrbayes.sourceforge.net/) or as Debian packages *mrbayes* and *mrbayes-mpi*; version used and recommended: 3.2.6)  
  - \[ future development: consider using [RevBayes](http://revbayes.github.io/) \]

- **MAD** for species tree rooting  
  ([R source code](https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip))

- **ALE/xODT** for gene tree / species tree reconciliation  
  (Install from [source code](https://github.com/ssolo/ALE); version used and recommended: 0.4; notably depends on [Bio++ libs](https://github.com/BioPP) (v2.2.0+))
  
### Required code libraries
- **R** (version 3, >=3.2.3 recommended) + packages:
  - ape
  - phytools
  - vegan
  - ade4
  - igraph
  - getopt
  - parallel
  - DBI, RSQLite
  - topGO (optional)
  - pvclust (optional)
  
- **Python** (version 2.7, >=2.7.13 recommended) + packages:
  - [sqlite3](https://docs.python.org/2/library/sqlite3.html) (standard package in Python 2.7)
  - [scipy/numpy](https://www.scipy.org/scipylib/download.html)
  - [tree2](https://github.com/flass/tree2)
  - [BioPython](http://biopython.org/wiki/Download)
  - [BCBio.GFF](https://pypi.org/project/bcbio-gff)
  - [Cython](https://pypi.org/project/Cython/)
  - [igraph](http://igraph.org/python/) (available as a Debian package)

### Other required software
- [sqlite3](https://www.sqlite.org) (available as a Debian package *sqlite3*)
- [LFTP](https://lftp.yar.ru/get.html) (available as a Debian package *lftp*)
- [(linux)brew](http://linuxbrew.sh/) (available as a Debian package *linuxbrew-wrapper*)
- [docker](https://www.docker.com/) (available as a Debian package *docker.io*)
- [JAVA Runtime (JDK 8.0)](https://openjdk.java.net) (available as Debian packages *openjdk-8-jdk* and *openjdk-8-jre*)
- [CD-HIT](https://cd-hit.org) (available as a Debian package *cd-hit*)
- [bioperl](https://bioperl.org) (available as a Debian package *bioperl*)
