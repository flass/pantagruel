FROM debian:stretch

RUN apt-get update && apt-get install -y git build-essential cmake gcc g++ \
        linuxbrew-wrapper lftp clustalo raxml libhmsbeagle1v5 mrbayes r-base-core \
        r-recommended r-cran-ape r-cran-ade4 r-cran-vegan r-cran-dbi r-cran-rsqlite \
        r-cran-igraph r-cran-getopt sqlite3 sqlite3-doc libmagick++-dev python \
        python-scipy python-numpy python-biopython python-biopython-sql python-igraph \
        cython bioperl mpi-default-bin mpi-default-dev mrbayes-mpi python-pip \
        openjdk-8-jdk openjdk-8-jre cd-hit

RUN echo 'source("https://bioconductor.org/biocLite.R") ; biocLite("topGO")' | R --vanilla
RUN apt-get install -y libxml2-dev libcurl4-openssl-dev
RUN echo 'install.packages(c("phytools","pvclust"), repos="https://pbil.univ-lyon1.fr/CRAN/")' | R --vanilla


RUN apt-get update && apt-get install -y libdatetime-perl libxml-simple-perl \
    libdigest-md5-perl default-jdk bioperl
RUN apt-get install -y wget

# Prokka install
RUN wget https://github.com/tseemann/prokka/archive/v1.12.tar.gz

RUN tar xzf v1.12.tar.gz && \
    rm v1.12.tar.gz && \
    cd /usr/bin && \
    ln -s /prokka-1.12/bin/prokka

RUN prokka --setupdb


#MMSEQ2
RUN apt-get install -y --allow-unauthenticated cmake xxd
RUN wget "https://github.com/soedinglab/MMseqs2/archive/1-c7a89.tar.gz" && \
    tar xvfz 1-c7a89.tar.gz
RUN cd MMseqs2-1-c7a89 && \
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

COPY scripts /pantagruel/scripts
COPY python_libs /pantagruel/python_libs
