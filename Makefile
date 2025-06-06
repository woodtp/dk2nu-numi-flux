CXX=g++
CFLAGS=-O3 -std=c++20 -Wall -Wextra -pedantic -c -fPIC -I${ROOT_INC}
DK2NU_BUILD=${PWD}/dk2nu/build

.PHONY: all dk2nu venv clean

all: lib dk2nu venv

lib: Weight.o
	$(CXX) -shared -o libWeight.so Weight.o

Weight.o: Weight.cc
	echo $(CFLAGS)
	$(CXX) $(CFLAGS) Weight.cc

dk2nu:
	echo "Building libdk2nuTree.so"
	git submodule update --init
	test -d $(DK2NU_BUILD) || mkdir -p $(DK2NU_BUILD)
	cmake -DWITH_GENIE=OFF -DWITH_TBB=OFF -B $(DK2NU_BUILD) -S $(DK2NU_BUILD)/.. && make -B -C $(DK2NU_BUILD) -j$(nproc)
	ln -sfv $(DK2NU_BUILD)/tree/{libdk2nuTree.rootmap,module.modulemap,libdk2nuTree_rdict.pcm} $(DK2NU_BUILD)/lib

venv:
	python -m venv venv
	. ./venv/bin/activate && pip install -r requirements.txt

clean:
	rm -fv libWeight.so Weight.o
	rm -rf ./venv
	rm -rf ./dk2nu/build/*
