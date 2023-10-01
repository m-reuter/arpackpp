#!/bin/bash
set -e

# The OpenBLAS version (git tag) to download.
version="0.3.24"

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

# Download the version specified above. In case the archive or extracted folder
# already exist, they will be re-used. Use the --cleanup option to remove the
# folder (and thus temporary build artifacts and CMake cache) after build.

if [ ! -f "OpenBLAS-${version}.tar.gz" ]; then
  wget -O OpenBLAS-${version}.tar.gz https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v${version}.tar.gz
fi

if [ ! -d "OpenBLAS-${version}" ]; then
  tar -xvf OpenBLAS-${version}.tar.gz > /dev/null
fi

cd OpenBLAS-${version}

# Set CMake install prefix options.
install_prefix_CONF=""
install_prefix_INST=""

if [ ! -z "$install_prefix" ]; then
  install_prefix_CONF="-D CMAKE_INSTALL_PREFIX=$install_prefix"
  install_prefix_INST="--prefix $install_prefix"
fi

cmake -B build -D BUILD_TESTING=OFF -D CMAKE_BUILD_TYPE=$build_type -D BUILD_SHARED_LIBS=$shared_libs $install_prefix_CONF
cmake --build build --config $build_type --parallel
cmake --install build $install_prefix_INST

cd ../../

if [ $cleanup -eq 1 ]; then
  echo "Cleaning up ..."
  rm -rf external/OpenBLAS-${version}/
fi
