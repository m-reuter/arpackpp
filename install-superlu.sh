#!/bin/sh
set -ex
sluversion="5.2.1"

if [ ! -d "external/SuperLU_${sluversion}" ]; then
  mkdir -p external
  cd external
  extdir=$PWD
  curl -o superlu_${sluversion}.tar.gz http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_${sluversion}.tar.gz
  tar -xvf superlu_${sluversion}.tar.gz > /dev/null
  rm superlu_${sluversion}.tar.gz

  # 5.0 needs some massaging:
  if [ "$sluversion" = "5.0" ]; then
    cd SuperLU_5.0
    cp make.inc make.inc.orig
    sed -i.bak '/PLAT/c\PLAT = ' make.inc
    sed -i.bak '/Dropbox/c\SuperLUroot \t= $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))' make.inc
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
  
  else # designed for current version 5.2.0:
  
    # rewrite required cmake version (it is unnecessary high)
    sed -i.bak 's/cmake_minimum_required(VERSION 2.8.12)/cmake_minimum_required(VERSION 2.8)/' ./SuperLU_${sluversion}/CMakeLists.txt
    # make use of out-of-source build with cmake
    mkdir SuperLU_${sluversion}-build
    cd SuperLU_${sluversion}-build
    # on my mac, I need to pass FC location:
    if [ "`uname`" = "Darwin" ] && [ -e "/opt/local/bin/gfortran-mp-4.9" ] ; then
      fcstr="-DCMAKE_Fortran_COMPILER=/opt/local/bin/gfortran-mp-4.9"
    fi
    # if we don't use system blas, use our local openblas 
    if [ ! "$BLAS" = "SYSTEM" ] ; then
      if [ ! -d "$extdir/OpenBLAS" ]; then
        echo 'ERROR: Local OpenBLAS does not exist, please install first'
        exit 1
       fi
       echo 'Using local openblas ..'
       blasstr="-DTPL_BLAS_LIBRARIES=$extdir/libopenblas.a"
    fi
    cmake $blasstr -Denable_blaslib=OFF -Denable_tests=OFF $fcstr ../SuperLU_${sluversion}
    make
    cd ../
    ln -s SuperLU_${sluversion}-build/SRC/libsuperlu.a ./
  fi
  
  ln -s SuperLU_${sluversion} SuperLU
  cd ../
  
else
  echo "external/SuperLU_${sluversion} already exists, using cached directory."
fi
# apt-get install libsuperlu-dev (still 4.3??)
