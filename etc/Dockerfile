FROM ubuntu:bionic

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y git build-essential cmake gcc g++ \
        linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes r-base-core \
        r-recommended r-cran-ape r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite \
        r-cran-igraph r-cran-getopt sqlite3 sqlite3-doc libmagick++-dev python \
        python-scipy python-numpy python-biopython python-biopython-sql python-igraph \
        cython mpi-default-bin mpi-default-dev mrbayes-mpi python-pip \
        openjdk-11-jdk openjdk-11-jre parallel libdw1 libdw-dev libgsl23 libgsl-dev

RUN echo 'source("https://bioconductor.org/biocLite.R") ; biocLite("topGO")' | R --vanilla
RUN apt-get install -y libxml2-dev libcurl4-openssl-dev
RUN echo 'install.packages(c("phytools","pvclust"), repos="https://pbil.univ-lyon1.fr/CRAN/")' | R --vanilla


# HOMEBREW
RUN apt-get install -y locales linuxbrew-wrapper
RUN localedef -i en_US -f UTF-8 en_US.UTF-8 \
   && useradd -m -s /bin/bash linuxbrew \
   && echo 'linuxbrew ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers
USER linuxbrew
WORKDIR /home/linuxbrew
ENV USER=linuxbrew \
    PATH=/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:$PATH \
    SHELL=/bin/bash
RUN HOMEBREW_NO_ANALYTICS=1 HOMEBREW_NO_AUTO_UPDATE=1 brew tap brewsci/bio
USER root
WORKDIR /

# PROKKA
RUN brew install brewsci/bio/prokka
RUN prokka --setupdb

# MMSEQS
RUN brew install mmseqs2

#PAL2NAL
RUN wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz && \
    tar -xzf pal2nal.v14.tar.gz && \
    chmod +x /pal2nal.v14/pal2nal.pl && \
    ln -s /pal2nal.v14/pal2nal.pl /usr/bin

##MAD
#RUN wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip && \
#    unzip ${SOFTWARE}/mad2-2.zip && \
#    chmod +x ${SOFTWARE}/mad/mad && \
#    ln -s ${SOFTWARE}/mad/mad /usr/bin

# PYTHON libs
RUN pip install bcbio-gff
RUN pip install bioscripts.convert

#INTERPROSCAN
#RUN (cd /usr/local && wget -O interproscan-5.35-74.0-64-bit.tar.gz ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.35-74.0/interproscan-5.35-74.0-64-bit.tar.gz)

# LSD
RUN wget -O /usr/bin/lsd https://github.com/tothuhien/lsd-0.3beta/releases/download/v0.3.3/lsd_unix && \
    chmod +x /usr/bin/lsd

RUN apt-get install -y rsync

# ALE
RUN apt-get install -y libboost-dev libboost-serialization-dev libboost-mpi-dev
RUN apt-get install -y libbpp-core-dev libbpp-phyl-dev libbpp-seq-dev libbpp-seq-omics-dev
RUN git clone https://github.com/ssolo/ALE
RUN cd ALE && mkdir build && cd build && cmake .. && make
RUN cp /ALE/build/bin/* /usr/bin

#
#COPY ./scripts /pantagruel/scripts
#COPY ./python_libs /pantagruel/python_libs
#COPY ./pantagruel /pantagruel
#COPY ./.git /pantagruel/.git
#ENV PATH="/pantagruel:/pantagruel/scripts:/pantagruel/scripts/pipeline:${PATH}"

