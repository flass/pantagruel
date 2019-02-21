FROM debian:stretch

COPY scripts /pantagruel/scripts
COPY python_libs /pantagruel/python_libs

RUN apt-get update && apt-get install -y git build-essential cmake gcc g++ linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes r-base-core r-recommended r-cran-ape r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite r-cran-igraph r-cran-getopt sqlite3 sqlite3-doc libmagick++-dev python python-scipy python-numpy python-biopython python-biopython-sql python-igraph cython bioperl mpi-default-bin mpi-default-dev mrbayes-mpi python-pip openjdk-8-jdk openjdk-8-jre cd-hit
