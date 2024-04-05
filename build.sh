source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

export DK2NU=$PWD/dk2nu
export DK2NU_INC=$DK2NU/tree
export DK2NU_LIB=$DK2NU/build/lib

setup root v6_28_10a -q e28:p3915:prof
setup cmake v3_27_4

BUILD_DIR=$DK2NU/build

mkdir $BUILD_DIR
cd $BUILD_DIR
cmake -DWITH_GENIE=OFF -DWITH_TBB=OFF ..
make -j$(nproc)

ln -sv $DK2NU/build/tree/{libdk2nuTree.rootmap,module.modulemap,libdk2nuTree_rdict.pcm} $DK2NU_LIB

cd ../..

python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
