# Installing Pantagruel and its dependencies

## Installing the Pantagruel scripts

Pantagruel pipeline is 'just' a bunch of scripts. What you simply need to do is to download them, that's about it. Here is how:

Under a Debian environment (e.g. Ubuntu), you can automatically install all dependencies using the script [install_dependencies.sh](https://github.com/flass/pantagruel/blob/master/install_dependencies.sh), following the indications below:  

Assuming you want to create a folder named `pantagruel_pipeline/` in the current working directory and install the whole software pipeline (the *pantagruel* package and its dependencies) in it, you should first do:
```sh
mkdir ./pantagruel_pipeline/
```
Then, you have to get the pantagruel scripts by dowloading the [archive of the last version on Github](https://github.com/flass/pantagruel/archive/master.zip) or use `git` to synchronize the repository (recomended for easier software update, especially during *Pantagruel* development phase!).
```sh
cd pantagruel_pipeline/
git clone --recurse-submodules https://github.com/flass/pantagruel.git
cd .. # return to parent folder
```

Finally, you need to install the dependencies. Because of course, the pipeline is more than just the master scripts, it's also all the fantastic software it realies on.

## Installing Dependencies 

### The container, worryless way :whale:

This is what is recommended, because that should always work. and that it provide a shared environment for us all in which you can report your bugs and I know what's going on.

First, you need Docker installed. On Debian-based systems, you should run:
```sh
sudo apt install docker.io
```

#### Downloading the docker image from Dockerhub
The simpler way to get the docker image for the pantagruel dependencies is to download it form the [Dockerhub](https://hub.docker.com/r/flass/pantagruel-dep) repository, where builds are made automatically with the latest code on the Gihub repository - magical!

Just run:
```sh
docker pull flass/pantagruel-dep:master-latest
```

#### Building the docker image
If for some reason, you want to build your own image rather than downloading it (you may have modified the source code for instance), do as follows:

Assuming you are still in the parent folder where you initially created the sub-folder `pantagruel_pipeline/`, run:

```sh
docker build -t panta pantagruel_pipeline/pantagruel
```

This should create a Docker image called `pantagruel-dep` that will be stored on your server. It contains all Pantagruel's dependecies (except InterProScan, see below) and the latest version of its source code. That's it!  

#### Using the docker image

Now, to run the pipeline, you will only need to create a docker container from that image. The pipeline scripts (and other executables) are available from the `$PATH` of the container, so you can use this syntax:  
```sh
docker run -u $UID:$UID -v $PWD:$PWD -w $PWD flass/pantagruel-dep:master-latest pantagruel
```
NB: if you built yourself the docker image, replace `flass/pantagruel-dep:master-latest` by `panta` or whatever else you used as value of option `docker build -t`.  

The `-v $PWD:$PWD` part of the command mounts your current directory `$PWD` as a volume onto the filesystem of the container, which otherwise is completely isolated from your host machine's own filesystem. The right bit of that part `:$PWD` indicates the mounting destination, which in that case is the same, meaning this location will have the same path in the container's filesystem, making it seemless to use the docker container vs. other modes of installing Pantagruel. Mounting a location of your own filesystem (it can be somewhere else than `$PWD`) thus allows you to write in it - so you can keep the output! The location chosen as root of your Pantgruel database (given as argument of `-r` option of the `init` task call) thus has to be included in that location. For instance, you can write the results in `/home/me/pantagruel_databases` if you mounted `/home/me` onto the container:  
```sh
docker run -u $UID:$UID -v /home/me:/home/me -w /home/me flass/pantagruel-dep:master-latest pantagruel -d nameofptgdb -r /home/me/pantagruel_databases init
```

If you want to use another version of Pantagruel source code (for instance because a bug fix or a new feature has been published), you don't have to rebuild the docker image! You can call an external version of pipeline through the container, provided the location of that code repository on your machine is included in the mounted volume:
```sh
docker run -u $UID:$UID -v $PWD:$PWD -w $PWD flass/pantagruel-dep:master-latest $PWD/pantagruel_pipeline/pantagruel/pantagruel
```

You can even alias this command so it's less ugly and you just need to call `pantagruel`:
```sh
alias pantagruel="docker run -u $UID:$UID -v $PWD:$PWD -w $PWD flass/pantagruel-dep:master-latest pantagruel_pipeline/pantagruel/pantagruel"
```

#### InterProScan/task 04 NOT included in docker image

Note that task `04` for InterProScan functional annotation is **NOT included** in the docker image, as InterProScan is bulky and frequent releases require regular manual re-installation.  
It can however be installed manually *in complement* of the docker image, and be called as an external program through the container; for this again you just need to make sure the executable `interproscan` is located somewhere within the folder mounted with `docker run` option `-v`.  
To this end, you can use the [install_interproscan.sh](https://github.com/flass/pantagruel/blob/master/install_interproscan.sh) script, using the same syntax as the `install_dependencies.sh` script (see above) but installing only InterProScan:  
```sh
pantagruel_pipeline/pantagruel/install_interproscan.sh pantagruel_pipeline/ $PWD/
```
This will download the last version of InterProScan, extract the (BIG!) Java library into `pantagruel_pipeline/`, and link the executable `interproscan.sh` to `$PWD/interproscan`; you can use any locations instead of `pantagruel_pipeline/` and `$PWD`, but remember they have to be located within the folder that will be mounted with `docker run -v`.

Then to indicate to Pantagruel where to find the executable, run the initial configuration command `pantagruel init` with the option `--path_to_interproscan $PWD/interproscan`.

### The scripted, fairly easy way

#### Using the `install_dependencies.sh` script 

After cloning the `pantagruel` code repository, you may run the installation script [install_dependencies.sh](https://github.com/flass/pantagruel/blob/master/install_dependencies.sh):  

Assuming you are still in the parent folder where you initially created the sub-folder `pantagruel_pipeline/`, run:
```sh
pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline/
```
This will download and build all dependencies in the `pantagruel_pipeline/` folder.

!!! Note that for this step, you (i.e. your linux user) need to have the sudo rights (be part of the sudo group, check with `grep sudo /etc/group`); however, **DO NOT** execute the installtion script with `sudo` (which would make all software owned by the root user).  

Optionally, you can also specify the folder where relevant executable files, including the main `pantagruel` executable, will be linked, or it will default to `~/bin/` (user-specific). The path to this folder will be **permatnently** added to your `$PATH` (via editing your ~/.bashrc).  
```sh
pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline/ pantagruel_pipeline/bin/
```

You may want not to install automatically all Debian packages (some could mess up with your local install), Brew and all its packages (same reason), Docker and all its packages (you may have a special deamon installed you don't want to be replaced) or InterProScan (it takes a lot of disk space). For this, you can use the options `--no-debian` `--no-brew` `--no-docker` and `--no-interpro`, respectively (anywhere after the first argument).  
```sh
pantagruel_pipeline/pantagruel/install_dependencies.sh pantagruel_pipeline/ --no-debian --no-interpro --no-brew --no-docker
```

#### Checking ALE functions

You may want to test that ALE commands work after installation (see [here](https://github.com/flass/pantagruel/blob/master/doc/installing_ALE.md#checking-it-works)), by typing the single commands `ALEml` and `ALEobserve` - preceeeded by `docker run --user $USER -v $PWD:$PWD -w $PWD` if using the docker image.



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
To compile ALE from source, Bio++ version >=2.2.0 is required. The version 2.4.1 can be found as a Debian package on Ubuntu 18.4 LTS (Bionic Beaver).  
```sh
sudo apt install libbpp-core-dev libbpp-phyl-dev libbpp-seq-dev libbpp-seq-omics-dev
```
Previous Ubuntu versions, such as 16.4 LTS (Xenial Xerius) have version 2.1.0, in which case Bio++ have to be compiled from source (much heavier, see [installing_ALE](https://github.com/flass/pantagruel/blob/master/doc/installing_ALE.md)).  
 
### Platform-independent software installation

in the following generic commands, the `${SOFTWARE}` and `${BINS}` environment variables would contain the path to folders where you want to install programs and where you want to link their executables, respectively.
For programs to be accessible for execution during by the pipeline, `${BINS}` should feature in your `${PATH}` variable definition:
```sh
export PATH=${PATH}:${BINS}
```
The command above need to be added to your `~/.bashrc` file for the programs to be permanently available.

#### Fetch Pantagruel pipeline and specific phylogenetic modules
```sh
# assuming you want to install this in a folder named ${SOFTWARE}:
cd ${SOFTWARE}/
git clone https://github.com/flass/pantagruel
cd pantagruel/ && git submodule init && git submodule update && cd -
```
You then need to add the Python modules to your `PYTHONPATH`:
```sh
export PYTHONPATH=${PYTHONPATH}:${SOFTWARE}/pantagruel/python_libs
```
The command above need to be added to your `~/.bashrc` file for the module to be permanently available.

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
ln -s ${SOFTWARE}/MAD/mad ${BINS}/
```
#### Fetch LSD program
```sh
wget https://github.com/tothuhien/lsd-0.3beta/releases/download/v0.3.3/lsd_unix
chmod +x ${SOFTWARE}/lsd_unix
ln -s ${SOFTWARE}/lsd_unix ${BINS}/lsd
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
  
- **LSD** for species tree rooting  
  ([source code](http://www.atgc-montpellier.fr/LSD/))
  
- **ALE/xODT** for gene tree / species tree reconciliation  
  (Install from [source code](https://github.com/ssolo/ALE); version used and recommended: 1.0; notably depends on [Bio++ libs](https://github.com/BioPP) (v2.2.0+))
  
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
- [bioperl](https://bioperl.org) (available as a Debian package *bioperl*)
