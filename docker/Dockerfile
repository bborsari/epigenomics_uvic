FROM debian:stable

LABEL maintainer "Diego Garrido Martín <diego.garrido@crg.eu>" \
      version "1.0" \
      description "Docker image for epigenomics course"

# Install basic tools
RUN apt-get update --fix-missing -qq && \
        apt-get install -y -q \
    bsdmainutils \
    bc \
    nano \ 
    wget \
    gawk \
    git-all \
    python \
    r-base \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir software

## Install BEDtools
RUN cd software \
    && wget --quiet https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz \
    && tar -zxvf bedtools-2.25.0.tar.gz \
    && cd bedtools2 \
    && make \
    && ln -s /software/bedtools2/bin/bedtools /usr/local/bin/bedtools

## Install bwtool
RUN cd software \
    && git clone https://github.com/CRG-Barcelona/libbeato.git \
    && git clone https://github.com/CRG-Barcelona/bwtool.git \
    && cd libbeato \
    && git checkout 0c30432af9c7e1e09ba065ad3b2bc042baa54dc2 \
       ## Latest commit is not working
    && ./configure \
    && make \
    && make install \
    && cd ../bwtool \
    && ./configure \
    && make \
    && make install

## Install bigBedToBed
RUN cd software \
    && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed \
    && chmod +x bigBedToBed \
    && ln -s /software/bigBedToBed /usr/local/bin/bigBedToBed

## Install R packages
RUN R -e 'install.packages( c("optparse", "ggplot2"), repos="http://cloud.r-project.org/");' 
