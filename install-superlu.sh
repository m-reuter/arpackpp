#!/bin/sh
set -ex
if [ ! -d "external/SuperLU" ]; then
  mkdir -p external
  cd external
  extdir=$PWD
  wget http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.0.tar.gz
  tar -xvf superlu_5.0.tar.gz
  cd SuperLU_5.0
  cp make.inc make.inc.orig
  sed -i '/PLAT/c\PLAT = ' make.inc
  sed -i '/Dropbox/c\SuperLUroot \t= $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))' make.inc
  if [ "$BLAS" = "SYSTEM" ] ; then
    echo 'Using system openblas'
  else
    if [ ! -d "$extdir/OpenBLAS" ]; then
      echo 'Local OpenBLAS does not exist, please install first'
    else
      echo 'Using local openblas'
      sed -i '/lib -lblas/c\BLASLIB \t= $(SuperLUroot)\/..\/OpenBLAS\/libopenblas.a  -lpthread' make.inc
    fi
  fi
  make lib
  cd ../
  ln -s SuperLU_5.0/lib/libsuperlu_5.0.a ./libsuperlu.a
  ln -s SuperLU_5.0 SuperLU
  cd ../
else
  echo 'external/SuperLU already exists, using cached directory.'
fi
# apt-get install libsuperlu-dev (still 4.3??)
