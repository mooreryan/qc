FROM mooreryan/build_essentials:0.1.0
LABEL maintainer="moorer@udel.edu"

ARG pipeline_version="0.7.0"
ARG software=/home/software

WORKDIR ${software}

# Install OpenJDK Java JRE.
RUN apt-get update && apt-get install -y openjdk-8-jre

# Bowtie2

## Download the Botiew2 binary.
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip
## Extract the archive.
RUN unzip bowtie2-2.3.5.1-linux-x86_64.zip
## Update the PATH variable with the Bowtie2 binary directory.
ENV PATH=${PATH}:${software}/bowtie2-2.3.5.1-linux-x86_64
## Remove the original archive.
RUN rm bowtie2-2.3.5.1-linux-x86_64.zip

# FixPairs

## Get the FixPairs source code.
RUN wget https://raw.githubusercontent.com/mooreryan/FixPairs/47064a2ca5709070c15df9098db2d97bfe937109/fix_pairs.cc
## It's a C++ program, so compile it.
RUN g++ -O2 -Wall --std=c++11 -o ${bin}/FixPairs fix_pairs.cc

# FLASH

## Download the flash binary.
RUN wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11-Linux-x86_64.tar.gz
## Extract the archive.
RUN tar xzf FLASH-1.2.11-Linux-x86_64.tar.gz
## Update the PATH variable with the flash binary directory.
ENV PATH=${PATH}:${software}/FLASH-1.2.11-Linux-x86_64
## Remove the original archive.
RUN rm FLASH-1.2.11-Linux-x86_64.tar.gz

# QC pipeline

## Download the source code
RUN wget https://github.com/mooreryan/qc/archive/v${pipeline_version}.tar.gz
## Extract the source code from the archive.
RUN tar xzf v${pipeline_version}.tar.gz
## Move into the directory of the source code.
WORKDIR qc-${pipeline_version}
## Change the permissions for the two pipeline scripts.
RUN chmod 755 qc.rb
RUN chmod 755 qc_multilib_wrapper.rb
## Install Ruby dependencies
RUN bundle install
## Move back to original working directory
WORKDIR ${software}
## Update the PATH variable with the QC pipeline directory.
ENV PATH=${PATH}:${software}/qc-${pipeline_version}
## Remove the original archive.
RUN rm v${pipeline_version}.tar.gz

# Set the working directory back to /home
WORKDIR /home
