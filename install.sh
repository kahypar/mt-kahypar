git submodule init
git submodule update

# Build KaHyPar Library Interface
cd external_tools/kahypar
git submodule init
git submodule update
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
sudo make install.library
cd ../../..

# Build Mt-KaHyPar
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make KaHyPar