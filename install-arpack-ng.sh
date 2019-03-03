#!/bin/sh
set -ex
mkdir -p external
cd external
git clone https://github.com/opencollab/arpack-ng.git
mkdir arpack-ng-build
cd arpack-ng-build
# on my mac, I need to pass FC location:
if [ "`uname`" = "Darwin" ] && [ -e "/opt/local/bin/gfortran-mp-4.9" ] 
then
  fcstr="-DCMAKE_Fortran_COMPILER=/opt/local/bin/gfortran-mp-4.9"
fi
cmake -D BLAS_goto2_LIBRARY=$PWD/../libopenblas.a $fcstr ../arpack-ng
make
cd ../
ln -s arpack-ng-build/libarpack.a ./
cd ../
