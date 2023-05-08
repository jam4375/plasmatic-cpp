# docker build -t plasmatic .

# Download base image ubuntu 22.04
FROM ubuntu:22.04

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

# Update Ubuntu Software repository
RUN apt update && apt upgrade -y && \
    apt install -y build-essential cmake git doxygen clang-tidy pkg-config wget python3 mpich libblas-dev liblapack-dev && \
    rm -rf /var/lib/apt/lists/* && apt clean

# Need to allow MPI running as root inside docker container
ENV OMPI_ALLOW_RUN_AS_ROOT 1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1
RUN mkdir /petsc && \
    cd /petsc && \
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.19.1.tar.gz && \
    tar xf petsc-3.19.1.tar.gz && \
    cd petsc-3.19.1 && \
    ./configure --prefix=/usr/local --with-fortran-bindings=0 \
                --with-debugging=no                          \
                --COPTFLAGS="-O3"                             \
                --CXXOPTFLAGS="-O3"                        && \
    make all && \
    make install

# Build the code
RUN mkdir /local
COPY . /local
ENV PKG_CONFIG_PATH /usr/local/lib/pkgconfig
RUN cd /local && rm -r build && mkdir -p build && cd build && cmake .. && make -j 16