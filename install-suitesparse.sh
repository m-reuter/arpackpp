#!/bin/bash
set -e

# The SuiteSparse version (git tag) to download.
version="7.6.1"

build_type="Release"
shared_libs="OFF"
cleanup=0

# Install into local folder "external", unless the --global-install option is used.
install_prefix="$(pwd)/external"

while [[ "$#" -gt 0 ]]; do
  case "${1:-}" in
    -g|--global-install)
      install_prefix=""
      shift 1
      ;;
    -s|--shared-libs)
      shared_libs="ON"
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

# Download the version specified above. In case the archive or extracted folder
# already exist, they will be re-used. Use the --cleanup option to remove the
# folder (and thus temporary build artifacts and CMake cache) after build.

source="https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${version}.tar.gz"
target="SuiteSparse-${version}.tar.gz"

if [ ! -f $target ]; then
  if [ -x "$(command -v curl)" ]; then
    curl -L -o $target $source
  elif [ -x "$(command -v wget)" ]; then
    wget -O $target $source
  else
    echo "Please install curl or wget!"
    exit 1
  fi
fi

if [ ! -d "SuiteSparse-${version}" ]; then
  tar -xvf $target > /dev/null
fi

cd SuiteSparse-${version}

# Set CMake install prefix options.
install_prefix_CONF=""
install_prefix_INST=""

if [ ! -z "$install_prefix" ]; then
  install_prefix_CONF="-D CMAKE_INSTALL_PREFIX=$install_prefix"
  install_prefix_INST="--prefix $install_prefix"
fi

cmake -B build -D SUITESPARSE_ENABLE_PROJECTS="cholmod;umfpack" -D CMAKE_BUILD_TYPE=$build_type -D BUILD_SHARED_LIBS=$shared_libs $install_prefix_CONF $local_BLAS
cmake --build build --config $build_type --parallel
cmake --install build $install_prefix_INST

cd ../../

if [ $cleanup -eq 1 ]; then
  echo "Cleaning up ..."
  rm -rf external/SuiteSparse-${version}/
fi
