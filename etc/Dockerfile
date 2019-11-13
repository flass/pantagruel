FROM ubuntu:bionic

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y git build-essential cmake gcc g++ \
        linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes r-base-core \
        r-recommended r-cran-ape r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite \
        r-cran-igraph r-cran-getopt sqlite3 sqlite3-doc libmagick++-dev python \
        python-scipy python-numpy python-biopython python-biopython-sql python-igraph \
        cython mpi-default-bin mpi-default-dev mrbayes-mpi python-pip \
        openjdk-8-jdk openjdk-8-jre cd-hit parallel

RUN echo 'source("https://bioconductor.org/biocLite.R") ; biocLite("topGO")' | R --vanilla
RUN apt-get install -y libxml2-dev libcurl4-openssl-dev
RUN echo 'install.packages(c("phytools","pvclust"), repos="https://pbil.univ-lyon1.fr/CRAN/")' | R --vanilla


#RUN brew update && brew doctor && brew install brewsci/bio/prokka


# RUN prokka --setupdb

RUN apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
RUN cpan Bio::Perl
RUN apt-get install ncbi-blast+
RUN git clone https://github.com/tseemann/prokka.git /prokka && cd /prokka && git checkout d6adb2c79ceedf # closest commit to 1.13 that fixes https://github.com/tseemann/prokka/issues/291
RUN cd /usr/bin && ln -s /prokka/bin/prokka
RUN prokka --setupdb


#MMSEQ2
RUN apt-get install -y --allow-unauthenticated cmake xxd
RUN wget "https://github.com/soedinglab/MMseqs2/archive/7-4e23d.tar.gz" && \
    tar xvfz 7-4e23d.tar.gz
RUN cd MMseqs2-7-4e23d && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. && \
    make && \
    make install && \
    cp bin/mmseqs /usr/bin

#PAL2NAL
RUN wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz && \
    tar -xzf pal2nal.v14.tar.gz && \
    chmod +x /pal2nal.v14/pal2nal.pl && \
    ln -s /pal2nal.v14/pal2nal.pl /usr/bin

#MAD
RUN wget https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip && \
    unzip ${SOFTWARE}/mad2-2.zip && \
    chmod +x ${SOFTWARE}/mad/mad && \
    ln -s ${SOFTWARE}/mad/mad /usr/bin

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

COPY ./scripts /pantagruel/scripts
COPY ./python_libs /pantagruel/python_libs
COPY ./pantagruel /pantagruel
COPY ./.git /pantagruel/.git
ENV PATH="/pantagruel:/pantagruel/scripts:/pantagruel/scripts/pipeline:${PATH}"

