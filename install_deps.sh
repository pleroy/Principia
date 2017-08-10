#!/bin/bash
set -euo pipefail

echo "Required prerequisites for build: build-essential clang libc++-dev libc++abi-dev monodevelop subversion git"
echo "Required runtime dependencies: libc++1"

#sudo apt-get install clang git unzip wget libc++-dev binutils make automake libtool curl cmake subversion

BASE_FLAGS="-fPIC -O3 -g"
# determine platform for bitness

PLATFORM=$(uname -s)
if [ "$PLATFORM" == "Darwin" ]; then
    C_FLAGS="$BASE_FLAGS -mmacosx-version-min=10.7 -arch x86_64"
elif [ "$PLATFORM" == "Linux" ]; then
	BITNESS=$(uname -m)
	if [ "$BITNESS" == "x86_64" ]; then
	  C_FLAGS="$BASE_FLAGS -m64"
        else
      	   C_FLAGS="$BASE_FLAGS -m32"
        fi
else
    C_FLAGS="$BASE_FLAGS"
fi

LD_FLAGS="$C_FLAGS -stdlib=libc++"
CXX_FLAGS="-std=c++14 $LD_FLAGS"

mkdir -p deps
cd deps

if [ ! -d "protobuf" ]; then
  git clone "https://github.com/mockingbirdnest/protobuf"
fi
pushd protobuf
git checkout master
git pull
./autogen.sh
if [ "$PLATFORM" == "Linux" ]; then
    ./autogen.sh # Really definitely needs to run twice on Ubuntu for some reason.
fi
./configure CC=clang CXX=clang++ CXXFLAGS="$CXX_FLAGS" LDFLAGS="$LD_FLAGS" LIBS="-lc++ -lc++abi"
make -j8
popd

if [ ! -d "glog" ]; then
  git clone "https://github.com/mockingbirdnest/glog"
fi
pushd glog
git checkout master
git pull
aclocal
automake
./configure CC=clang CXX=clang++ CFLAGS="$C_FLAGS" CXXFLAGS="$CXX_FLAGS" LDFLAGS="$LD_FLAGS" LIBS="-lc++ -lc++abi"
make -j8
popd

# googlemock/googletest don't need to be compiled
if [ ! -d "googlemock" ]; then
  git clone "https://github.com/mockingbirdnest/googlemock"
fi
pushd googlemock
git checkout master
git pull
popd

if [ ! -d "googletest" ]; then
  git clone "https://github.com/mockingbirdnest/googletest"
fi
pushd googletest
git checkout master
git pull
popd

if [ ! -d "eggsperimental_filesystem" ]; then
  git clone "https://github.com/Norgg/eggsperimental_filesystem.git"
fi
pushd eggsperimental_filesystem
git checkout master
git pull
popd

# Optional doesn't need to be compiled either
if [ ! -d "Optional" ]; then
  git clone "https://github.com/mockingbirdnest/Optional.git"
fi
pushd Optional
git checkout master
git pull
popd

# TODO(egg): This probably needs to be compiled
if [ ! -d "benchmark" ]; then
  git clone "https://github.com/mockingbirdnest/benchmark"
fi
pushd benchmark
git checkout master
git pull
cmake -DCMAKE_C_COMPILER:FILEPATH=`which clang` -DCMAKE_CXX_COMPILER:FILEPATH=`which clang++` -DCMAKE_C_FLAGS="${C_FLAGS}" -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DCMAKE_LD_FLAGS="${LD_FLAGS}"
make -j8
popd
