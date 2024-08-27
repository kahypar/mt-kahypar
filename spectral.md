# julia

This branch uses [julia](https://julialang.org). It is tested under version 1.10.4.

The required julia packages are listed in `mt-kahypar/partition/refinement/spectral/solvers/julia/dependencies.txt`.

# K-SpecPart dependencies

To use spectral refinement, you'll need the following steps to install the dependencies to a local folder (just follow the respective official instructions for installation with root):

## spdlog

 - `git clone --branch v1.8.1 https://github.com/gabime/spdlog.git spdlog-1.8.1`
 - install via `cmake` with install prefix

## lemon

 - download from https://lemon.cs.elte.hu/trac/lemon
 - install via `cmake` with install prefix

## or-tools

 - `git clone https://github.com/google/or-tools`
 - install via `cmake` with install prefix and `BUILD_DEPS=ON`

## boost

 - download from https://www.boost.org/users/history/version_1_80_0.html
 - install with `bootstrap.sh --prefix=/install/prefix/`and `b2 install`

## tcl

 - download from http://prdownloads.sourceforge.net/tcl/tcl8.6.14-src.tar.gz
 - installation:
```
cd unix
configure --prefix=/install/prefix/
make install install-private-headers
cd /install/prefix/bin
ln -s tclsh8.6 tclsh
chmod +x tclsh
```

## swig

 - download from https://swig.org/svn.html
 - install via `autogen.sh`, `configure` (with install prefix) and `make`

## bison

 - download from http://ftp.gnu.org/gnu/bison/bison-3.8.tar.gz
 - install via `configure` with install prefix and `make`

## flex

 - download from https://github.com/westes/flex/files/981163/flex-2.6.4.tar.gz
 - install via `configure` with install prefix and `make`

## openroad

 - download from https://github.com/juliannz/TritonPart_OpenROAD.git
 - installation:
```
cmake .. -DENABLE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=/install/prefix/ -DTCL_LIBRARY=/install/prefix/lib/libtcl8.6.so -DUSE_SYSTEM_OR_TOOLS=ON
rm /install/prefix/bin/openroad
ln build/src/openroad /install/prefix/bin/openroad
```

## GKlib

 - download from https://github.com/KarypisLab/GKlib
 - install via `make` with install prefix

## METIS
 - download from https://github.com/KarypisLab/METIS
 - install via `make` with install prefix

## hmetis

 - download from http://glaros.dtc.umn.edu/gkhome/fetch/sw/hmetis/hmetis-1.5-linux.tar.gz
 - installation:
```
cp hmetis khmetis shmetis /install/prefix/bin
cp libhmetis.a /install/prefix/lib
cd /install/prefix/bin
chmod +x hmetis khmetis shmetis
```

## ilp_partitioner

 - `git clone https://github.com/juliannz/ilp_partitioner.git`
 - `cmake .. -DILP_USE_CPLEX=OFF -DCMAKE_INSTALL_PREFIX=/install/prefix/`