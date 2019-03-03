#!/bin/sh
#set -ex
set -x
mkdir -p external
cd external
rm -rf OpenBLAS
git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS
# on my mac, I need to pass FC location:
if [ "`uname`" = "Darwin" ] && [ -e "/opt/local/bin/gfortran-mp-4.9" ] 
then
  fcstr="FC=/opt/local/bin/gfortran-mp-4.9"
fi
make $fcstr
if [ $? -eq 0 ] && [ -e "libopenblas.a" ]
then
  echo "libopenblas.a successfully created file"
else
  # could be new CPU, but old binutils (e.g. Centos6)
  echo "Could not make lib, trying with NO_AVX2=1"
  make clean
  make NO_AVX2=1 $fcstr
  if [ $? -eq 0 ] && [ -e "libopenblas.a" ]
  then
    echo "libopenblas.a successfully created file"
  else
    echo "Could not make lib" >&2
    make clean
    cd ../../
    exit $?
  fi
fi
cd ../
rm -f libopenblas.a libopenblas.so
ln -s OpenBLAS/libopenblas.a ./
#ln -s OpenBLAS/libopenblas.so ./
cd ../
