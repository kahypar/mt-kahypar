#!/bin/bash
source scripts/submodule_heads.sh
ROOT=${PWD}

# Initialize GOOGLETEST
[ ! "$(ls -A external_tools/googletest)" ] &&
git clone https://github.com/google/googletest.git external_tools/googletest &&
cd external_tools/googletest && git checkout ${GOOGLETEST_HEAD} && cd ${ROOT}

# Initialize KAHYPAR
[ ! "$(ls -A external_tools/kahypar)" ] &&
git clone https://github.com/SebastianSchlag/kahypar.git external_tools/kahypar &&
cd external_tools/kahypar && git checkout ${KAHYPAR_HEAD} && cd ${ROOT}

# Initialize WHFC
[ ! "$(ls -A external_tools/WHFC)" ] &&
git clone https://github.com/larsgottesbueren/WHFC.git external_tools/WHFC &&
cd external_tools/WHFC && git checkout ${WHFC_HEAD} && cd ${ROOT}

# Initialize PYBIND11
[ ! "$(ls -A python/pybind11)" ] &&
git clone https://github.com/pybind/pybind11.git python/pybind11 &&
cd python/pybind11 && git checkout ${PYBIND11_HEAD} && cd ${ROOT}
