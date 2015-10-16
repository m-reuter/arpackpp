#!/bin/sh
set -ex
mkdir -p external
cd external
extdir=$PWD
wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.5.tar.gz
tar -xvf SuiteSparse-4.4.5.tar.gz
cd SuiteSparse
# get metis
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
tar -xvf metis-4.0.3.tar.gz
mv metis-4.0.3 metis-4.0
cd metis-4.0
cp Makefile.in Makefile.in.orig
sed -i '/CC = cc/c\CC = gcc' Makefile.in
sed -i '/OPTFLAGS = -O2/c\OPTFLAGS = -O3' Makefile.in
# if we use local static openblas:
cd ../SuiteSparse_config
cp SuiteSparse_config.mk SuiteSparse_config.mk.orig
sed -i "/BLAS = -lopenblas/c\BLAS = $extdir/libopenblas.a -lpthread" SuiteSparse_config.mk
cd ..
# build SuiteSparse (incl. metis)
make
cd ../../
