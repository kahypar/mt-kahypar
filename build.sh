#!/bin/bash

function get_num_cores {
  if [[ $(uname) == "Linux" ]]; then grep -c ^processor /proc/cpuinfo; fi
  if [[ $(uname) == "Darwin" ]]; then sysctl -n hw.ncpu; fi
}

ROOT=${PWD}
if [ -d .git ]; then
  # Mt-KaHyPar is build from a git repository
  git submodule update --init;
else
  # Mt-KaHyPar is build from a release archive
  # which does not include submodules
  ./scripts/checkout_submodules.sh
fi

CMAKE_BINARY="/home/tobias/cmake-3.20.3-linux-x86_64/bin/cmake"
CMAKE_COMMANDS=$1
if [ ! -f build/Makefile ]; then
  mkdir -p build
fi

cd build && $CMAKE_BINARY .. -DCMAKE_BUILD_TYPE=Release $CMAKE_COMMANDS && cd ${ROOT}
$CMAKE_BINARY  --build build --parallel "$(get_num_cores)" --target MtKaHyPar

