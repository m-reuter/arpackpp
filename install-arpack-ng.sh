#!/bin/bash
set -e

install_prefix="--prefix $(pwd)/external"
cleanup=0

while [[ "$#" -gt 0 ]]; do
  case "${1:-}" in
    -g|--global-install)
      install_prefix=""
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

# In case the target directory already exists, try to pull the
# latest changes. Otherwise, clone the repository.

if [ -d "arpack-ng" ] && [ -n "$(ls -A arpack-ng)" ]; then
  cd arpack-ng
  git pull
else
  git clone https://github.com/opencollab/arpack-ng.git
  cd arpack-ng
fi

cmake -B build -D TESTS=OFF
cmake --build build --config Release --parallel
cmake --install build $install_prefix

cd ../../

if [ $cleanup -eq 1 ]; then
  echo "Cleaning up ..."
  rm -rf external/arpack-ng/
fi
