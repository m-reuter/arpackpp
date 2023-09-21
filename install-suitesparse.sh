#!/bin/bash
set -e

cleanup=0

while [[ "$#" -gt 0 ]]; do
  case "${1:-}" in
    -l|--local-install)
      install_prefix="--prefix $(pwd)/external"
      echo "Local install enabled"
      shift 1
      ;;
    -c|--cleanup)
      cleanup=1
      shift 1
      ;;
  esac
done

# In case you want CMake to use a different compiler than the system default, use
#
#    export CC=opt/local/bin/your-c-compiler
#    export FC=opt/local/bin/your-fortran-compiler

mkdir -p external
cd external

if [ -d "SuiteSparse" ] && [ -n "$(ls -A SuiteSparse)" ]; then
  cd SuiteSparse
  git pull
else
  git clone https://github.com/DrTimothyAldenDavis/SuiteSparse.git
  cd SuiteSparse
fi

for target in (SuiteSparse_config AMD CAMD COLAMD CCOLAMD CHOLMOD UMFPACK); do
  cd $target
  cmake -B build
  cmake --build build --parallel
  cmake --install build $install_prefix
  cd ../
done

cd ../../

if [ $cleanup -eq 1 ]; then
  echo "Cleaning up ..."
  rm -rf external/SuiteSparse/
fi
