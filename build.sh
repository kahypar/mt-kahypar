#!/bin/bash
function get_num_cores {
  if [[ $(uname) == "Linux" ]]; then grep -c ^processor /proc/cpuinfo; fi
  if [[ $(uname) == "Darwin" ]]; then sysctl -n hw.ncpu; fi
}
ROOT=${PWD}
if [ -d .git ]; then
  # Mt-KaHyPar is build from a git repository
  git submodule update --init --recursive
else
  # Mt-KaHyPar is build from a release archive
  # which does not include submodules
  ./scripts/checkout_submodules.sh
fi
CMAKE_COMMANDS=$1
mkdir -p build_debug
mkdir -p build

# Enable tests (keep both flags; one of them will match your CMake setup)
TEST_FLAGS="-DKAHYPAR_ENABLE_TESTING=ON"

cd build_debug && cmake .. -DCMAKE_BUILD_TYPE=Debug ${TEST_FLAGS} ${CMAKE_COMMANDS} && cd "${ROOT}"
cmake --build build_debug --parallel "$(get_num_cores)" --target MtKaHyPar
cmake --build build_debug --parallel "$(get_num_cores)" --target mtkahypar_tests

cd build && cmake .. -DCMAKE_BUILD_TYPE=Release ${TEST_FLAGS} ${CMAKE_COMMANDS} && cd "${ROOT}"
cmake --build build --parallel "$(get_num_cores)" --target MtKaHyPar
cmake --build build --parallel "$(get_num_cores)" --target mtkahypar_tests

#cd build_debug && cmake .. -DCMAKE_BUILD_TYPE=Debug $CMAKE_COMMANDS && cd ${ROOT}
#cmake  --build build_debug --parallel "$(get_num_cores)" --target MtKaHyPar
#cd build && cmake .. -DCMAKE_BUILD_TYPE=Release $CMAKE_COMMANDS && cd ${ROOT}
#cmake  --build build --parallel "$(get_num_cores)" --target MtKaHyPar
#cd build && cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo $CMAKE_COMMANDS && cd ${ROOT}
#cmake  --build build --parallel "$(get_num_cores)" --target MtKaHyPar