FROM debian:stretch as intermediate

# install git
RUN apt-get update
RUN apt-get install -y git

# add credentials on build
ARG SSH_PRIVATE_KEY
RUN mkdir /root/.ssh/
RUN echo "${SSH_PRIVATE_KEY}" > /root/.ssh/id_rsa
RUN chmod 600 /root/.ssh/id_rsa

# make sure your domain is accepted
RUN touch /root/.ssh/known_hosts
RUN ssh-keyscan github.com >> /root/.ssh/known_hosts

RUN git clone git@github.com:flass/pantagruel.git


FROM pveber/bistro-base:stretch

COPY --from=intermediate /pantagruel /pantagruel
RUN chmod a+x /pantagruel/scripts/* && cp /pantagruel/scripts/* /usr/bin
RUN apt-get update
RUN apt-get install -y python-scipy r-base python-biopython sqlite3 git
RUN git clone https://github.com/flass/tree2.git
ENV PYTHONPATH /tree2

RUN wget "http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz" && \
    tar xvfz pal2nal.v14.tar.gz && \
    cd pal2nal.v14 && \
    cp pal2nal.pl /usr/bin
