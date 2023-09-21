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

if [ -d "superlu" ] && [ -n "$(ls -A superlu)" ]; then
  cd superlu
  git pull
else
  git clone https://github.com/xiaoyeli/superlu.git
  cd superlu
fi

cmake -B build
cmake --build build --parallel
cmake --install build $install_prefix

cd ../../

if [ $cleanup -eq 1 ]; then
  echo "Cleaning up ..."
  rm -rf external/superlu/
fi
