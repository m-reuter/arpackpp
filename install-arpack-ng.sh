#!/bin/sh
set -ex
mkdir -p external
cd external
git clone https://github.com/opencollab/arpack-ng.git
mkdir arpack-ng-build
cd arpack-ng-build
cmake -D BLAS_goto2_LIBRARY=../external/libopenblas.a ../arpack-ng
make
cd ../
ln -s arpack-ng-build/libarpack.a ./
cd ../
