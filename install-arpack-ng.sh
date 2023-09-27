#!/bin/bash
set -e

build_type="-D BUILD_SHARED_LIBS=OFF"
install_prefix_CONF="-D CMAKE_INSTALL_PREFIX=$(pwd)/external"
install_prefix_INST="--prefix $(pwd)/external"
cleanup=0

while [[ "$#" -gt 0 ]]; do
  case "${1:-}" in
    -g|--global-install)
      install_prefix_CONF=""
      install_prefix_INST=""
      echo "Global install (needs root access)"
      shift 1
      ;;
    -c|--cleanup)
      cleanup=1
      shift 1
      ;;
  esac
done

mkdir -p external
cd external

# Local BLAS lib (ignored if not exist)
lib_BLAS="$(pwd)/lib/libopenblas.a"

local_BLAS=""

if [[ -f "$lib_BLAS" ]]; then
  local_BLAS="-D BLAS_LIBRARIES=${lib_BLAS} -D LAPACK_LIBRARIES=${lib_BLAS}"
fi

# In case the target directory already exists, try to pull the
# latest changes. Otherwise, clone the repository.

if [ -d "arpack-ng" ] && [ -n "$(ls -A arpack-ng)" ]; then
  cd arpack-ng
  git pull
else
  git clone https://github.com/opencollab/arpack-ng.git
  cd arpack-ng
fi

cmake -B build -D TESTS=OFF $build_type $install_prefix_CONF $local_BLAS
cmake --build build --config Release --parallel
cmake --install build $install_prefix_INST

cd ../../

if [ $cleanup -eq 1 ]; then
  echo "Cleaning up ..."
  rm -rf external/arpack-ng/
fi
