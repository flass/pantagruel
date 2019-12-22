## Using Docker to install ALE
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
## Installing ALE from source code

This is more cumbersome, but safer when installing this software on a computer that has other functions, especially on  a large production server, or to benefit from the last development version. See [ALE's own installation page](https://github.com/ssolo/ALE/blob/master/INSTALL.md) for more details and updates.

### Getting the required libs with Debian Packages 

#### Boost libraries
```sh
sudo apt install libboost-dev libboost-serialization-dev libboost-mpi-dev
```
#### Bio++ libraries version >=2.2.0 are required
The version 2.4.1 can be found as a Debian package on Ubuntu 18.4 LTS (Bionic Beaver). Previous Ubuntu versions, such as 16.4 LTS (Xenial Xerius) have version 2.1.0, in which case Bio++ have to be compiled from source (much heavier).  
Using Debian packages:  
```sh
sudo apt install libbpp-core-dev libbpp-phyl-dev libbpp-seq-dev libbpp-seq-omics-dev
```

### Building the required libs from source

I won't detail building Boost C++ libs from source as there is no reason to do that; however, would you want to do so, please refer to the [Boost website](https://www.boost.org/).

#### Compiling Bio++ libraries from source

As explained above, this may be required as not all Debian system come with the latest and required version. Note however that compiling these C++ libraries is computationally heavy and long (over an hour).
  
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

### Compiling ALE from source

Finally install ALE:  
```sh
git clone https://github.com/ssolo/ALE
mkdir build
cd build
cmake ..
make
```

## Checking it works!

To verify that the ALE commands are available in your system (i.e. that the location of excutable is listed in your `$PATH`) and that they function correctly, you should simply type the commands in your shell, and get these usage messages back:  

```sh
$ ALEml
ALEml using ALE v0.4
usage:
 ./ALEml species_tree.newick gene_tree_sample.ale  [samples] [gene_name_separator]
$ ALEobserve 
ALEobserve using ALE v0.4
usage:
 ./ALEobserve gene_tree_sample.newicks [burnin=0]
```
