export DK2NU=$PWD/dk2nu
export DK2NU_INC=$DK2NU/tree
export DK2NU_LIB=$DK2NU/build/lib

source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup root v6_28_10a -q e28:p3915:prof
setup cmake v3_27_4
[[ -d ./venv ]] && source ./venv/bin/activate
