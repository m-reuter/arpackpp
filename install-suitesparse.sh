#!/bin/bash
set -ex
mkdir -p external
cd external
extdir=$PWD
curl -O http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.5.tar.gz
tar -xvf SuiteSparse-4.4.5.tar.gz > /dev/null
rm SuiteSparse-4.4.5.tar.gz
cd SuiteSparse
# get metis
curl -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
tar -xvf metis-4.0.3.tar.gz > /dev/null
rm metis-4.0.3.tar.gz
mv metis-4.0.3 metis-4.0
cd metis-4.0
cp Makefile.in Makefile.in.orig
if [[ "$OSTYPE" == "darwin"* ]]; then
  sed -i "" "s/CC = cc/CC = gcc/g" Makefile.in
  sed -i "" "s/OPTFLAGS = -O2/OPTFLAGS = -O3/g" Makefile.in
else
  sed -i '/CC = cc/c\CC = gcc' Makefile.in
  sed -i '/OPTFLAGS = -O2/c\OPTFLAGS = -O3' Makefile.in
fi
# if we use local static openblas:
cd ../SuiteSparse_config
cp SuiteSparse_config.mk SuiteSparse_config.mk.orig
if [[ "$OSTYPE" == "darwin"* ]]; then
  cp SuiteSparse_config_Mac.mk SuiteSparse_config.mk
else
  if [ "$BLAS" = "SYSTEM" ] ; then
    echo 'Using system openblas'
  else
    if [ ! -d "$extdir/OpenBLAS" ]; then
      echo 'Local OpenBLAS does not exist, please install first'
    else
      echo 'Using local openblas'
      sed -i "/BLAS = -lopenblas/c\BLAS = ${extdir}\/OpenBLAS\/libopenblas.a -lpthread" SuiteSparse_config.mk
    fi
  fi
fi
cd ..
# build SuiteSparse (incl. metis)
make library
cd ../../
