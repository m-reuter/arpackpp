#!/bin/sh
set -ex
mkdir -p external
cd external
wget http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_5.0.tar.gz
tar -xvf superlu_5.0.tar.gz
cd SuperLU_5.0
cp make.inc make.inc.orig
sed -i '/PLAT/c\PLAT = ' make.inc
sed -i '/Dropbox/c\SuperLUroot \t= $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))' make.inc
sed -i '/lib -lblas/c\BLASLIB \t= $(SuperLUroot)\/..\/OpenBLAS\/libopenblas.a  -lpthread' make.inc
make
cd ../
ln -s SuperLU_5.0/lib/libsuperlu_5.0.a ./libsuperlu.a
ln -s SuperLU_5.0 SuperLU
cd ../
# apt-get install libsuperlu-dev (still 4.3??)
