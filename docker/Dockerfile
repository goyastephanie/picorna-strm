FROM ubuntu:20.04

ARG SAMTOOLS_VERSION=1.21 
ARG BCFTOOLS_VERSION=1.15
ARG HTSLIB_VERSION=1.15
ARG TRIMMOMATIC_VERSION=0.39
ARG PICARD_VERSION=2.26.11
ARG SEQKIT_VERSION=2.1.0
ARG IVAR_VERSION=1.3.1
ARG SEQTK_VERSION=1.3
ARG PANDEPTH_VERSION=2.25

ARG BWA_MEM_VERSION=2.2.1
ARG BWA_VERSION=0.7.19

ENV DEBIAN_FRONTEND=noninteractive 
ENV PATH=$PATH:/opt/bbmap
ENV TRIMMOMATIC=/usr/local/bin/trimmomatic.jar
ENV PICARD=/usr/local/bin/picard.jar
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

RUN apt-get update && apt-get install -y \
    openjdk-11-jdk bc python3 wget libncurses5-dev bzip2 gcc g++ make zlib1g-dev libbz2-dev liblzma-dev unzip autoconf

# ---------- 
# Samtools
RUN cd /tmp \
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install
# ---------- 

# Trimmomatic
RUN cd /tmp \
    && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
    && unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
    && mv Trimmomatic-${TRIMMOMATIC_VERSION}/trimmomatic-${TRIMMOMATIC_VERSION}.jar /usr/local/bin/trimmomatic.jar
#---------- 

# BWA-MEM2
RUN cd /tmp \
    && wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM_VERSION}/bwa-mem2-${BWA_MEM_VERSION}_x64-linux.tar.bz2 \
    && tar jxf bwa-mem2-${BWA_MEM_VERSION}_x64-linux.tar.bz2 \
    && mv bwa-mem2-2.2.1_x64-linux/* /usr/local/bin/
#---------- 

# BWA 
RUN cd /tmp \
    && wget https://github.com/lh3/bwa/archive/refs/tags/v${BWA_VERSION}.tar.gz \
    && tar -xvzf v${BWA_VERSION}.tar.gz \
    && ls \
    && cd bwa-${BWA_VERSION} && make \
    && mv bwa /usr/local/bin/
#---------- 

# iVar and HTSLib
RUN cd /tmp \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar xjvf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure \
    && make \
    && make install \
    && cd /tmp \
    && wget https://github.com/andersen-lab/ivar/archive/refs/tags/v${IVAR_VERSION}.tar.gz \
    && tar xzvf v${IVAR_VERSION}.tar.gz \
    && cd ivar-${IVAR_VERSION} \
    && ./autogen.sh \
    && ./configure \
    && make \
    && make install
#----------

#PanDepth
RUN cd /tmp \
    && wget https://github.com/HuiyangYu/PanDepth/releases/download/v${PANDEPTH_VERSION}/PanDepth-${PANDEPTH_VERSION}-Linux-x86_64.tar.gz \
    && tar zvxf PanDepth-${PANDEPTH_VERSION}-Linux-x86_64.tar.gz \
    && cd PanDepth-${PANDEPTH_VERSION}-Linux-x86_64 \
    && mv pandepth /usr/local/bin
# ----------

# Clean up
RUN apt autoremove -y wget libncurses5-dev bzip2 gcc g++ make zlib1g-dev libbz2-dev liblzma-dev unzip autoconf \
    && rm -rf /var/lib/apt/lists/* /tmp/*

CMD ["bash"]
