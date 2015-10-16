#!/bin/sh
set -ex
mkdir -p external
cd external
git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS
make
cd ../
ln -s OpenBLAS/libopenblas.a ./
#ln -s OpenBLAS/libopenblas.so ./
cd ../
