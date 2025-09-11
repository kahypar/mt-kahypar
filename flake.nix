{
  description = "Dev environment for Mt-KaHyPar - Multi-Threaded Karlsruhe Graph and Hypergraph Partitioner (CLion + CMake)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        stdenv = pkgs.stdenv;
      in
      {
        devShells.default = pkgs.mkShell {
          name = "mt-kahypar-dev";

          nativeBuildInputs = with pkgs; [
            cmake
            ninja
            pkg-config
            gdb
            ccache
            clang-tools 
            boost
            tbb_2022_0
          ];

          buildInputs = with pkgs; [
            boost
            tbb_2022_0
            zlib
            fmt
            hwloc
            python3
            lld
          ];

          NIX_ENFORCE_NO_NATIVE=false;
          
          shellHook = ''
            export CC=${stdenv.cc}/bin/cc
            export CXX=${stdenv.cc}/bin/c++
            export CMAKE_PREFIX_PATH=${pkgs.tbb_2022_0.dev}:${pkgs.boost.dev}:$CMAKE_PREFIX_PATH
            echo "mt-kahypar dev shell ready."
          '';
        };
      });
}
