#!/bin/bash

# git submodule update --init
#
# PROJECT_ROOT=$PWD
# DK2NU=$PROJECT_ROOT/dk2nu
# DK2NU_INC=$DK2NU/tree
# DK2NU_LIB=$DK2NU/build/lib
#
# BUILD_DIR=$DK2NU/build

# if [ ! -f $DK2NU_LIB/libdk2nuTree.so ]; then
#     echo "Building libdk2nuTree.so"
#
#     mkdir $BUILD_DIR
#     cd $BUILD_DIR
#     cmake -DWITH_GENIE=OFF -DWITH_TBB=OFF ..
#     make -j$(nproc)
#
#     ln -sv $DK2NU/build/tree/{libdk2nuTree.rootmap,module.modulemap,libdk2nuTree_rdict.pcm} $DK2NU_LIB
#     cd $PROJECT_ROOT
# else
#     echo "libdk2nuTree.so found! ($DK2NU_LIB/libdk2nuTree.so)"
# fi

LIBWEIGHT="libWeight.so"

LIBRARY_MTIME=$(stat -c %Y $LIBWEIGHT)
WEIGHT_MTIME=$(stat -c %Y ./Weight.cc)

if [ $LIBRARY_MTIME -lt $WEIGHT_MTIME ]; then
    echo "Rebuilding $LIBWEIGHT"
    rm -f $LIBWEIGHT
fi

CXX=g++
CFLAGS="-O3 -Wall -Wextra -pedantic -c -fPIC $(root-config --cflags)"
$CXX $CFLAGS Weight.cc
$CXX -shared -o libWeight.so Weight.o

[[ ! -d "./venv" ]] && python -m venv ./venv
source ./venv/bin/activate
pip install -U pip
pip install -r requirements.txt
