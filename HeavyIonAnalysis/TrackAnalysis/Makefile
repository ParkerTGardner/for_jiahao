ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

MKDIR_BIN=mkdir -p $(PWD)/bin

all: mkdirBin  bin/PYTHIA_gen.exe bin/SHERPA_gen.exe 

mkdirBin:
	$(MKDIR_BIN)

bin/PYTHIA_gen.exe: src/PYTHIA_gen.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen.exe src/PYTHIA_gen.C
bin/SHERPA_gen.exe: src/SHERPA_gen.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/SHERPA_gen.exe src/SHERPA_gen.C


clean:
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/src/*.so
	rm -f $(PWD)/src/*.d
	rm -f $(PWD)/src/*.pcm
	rm -f $(PWD)/bin/*.exe
	rmdir bin
