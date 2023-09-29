ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

MKDIR_BIN=mkdir -p $(PWD)/bin

all: mkdirBin   bin/PYTHIA_gen_data_di_2.exe bin/PYTHIA_gen_data_di_3.exe bin/PYTHIA_gen_di_2.exe bin/PYTHIA_gen_di_3.exe 
#  bin/PYTHIA_gen_di.exe  bin/PYTHIA_gen_di_copy_2.exe bin/PYTHIA_gen_di_copy_3.exe  bin/PYTHIA_gen_di_copy.exe 

 
mkdirBin:
	$(MKDIR_BIN)

bin/PYTHIA_gen_data_di_2.exe: src/PYTHIA_gen_data_di_2.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen_data_di_2.exe src/PYTHIA_gen_data_di_2.C
bin/PYTHIA_gen_data_di_3.exe: src/PYTHIA_gen_data_di_3.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen_data_di_3.exe src/PYTHIA_gen_data_di_3.C
#bin/PYTHIA_gen_di.exe: src/PYTHIA_gen_di.C
#	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen_di.exe src/PYTHIA_gen_di.C
bin/PYTHIA_gen_di_2.exe: src/PYTHIA_gen_di_2.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen_di_2.exe src/PYTHIA_gen_di_2.C
bin/PYTHIA_gen_di_3.exe: src/PYTHIA_gen_di_3.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen_di_3.exe src/PYTHIA_gen_di_3.C
#bin/PYTHIA_gen_di_copy.exe: src/PYTHIA_gen_di_copy.C
#	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen_di_copy.exe src/PYTHIA_gen_di_copy.C
#bin/PYTHIA_gen_di_copy_2.exe: src/PYTHIA_gen_di_copy_2.C
#	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen_di_copy_2.exe src/PYTHIA_gen_di_copy_2.C
#bin/PYTHIA_gen_di_copy_3.exe: src/PYTHIA_gen_di_copy_3.C
#	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/PYTHIA_gen_di_copy_3.exe src/PYTHIA_gen_di_copy_3.C
clean:
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/src/*.so
	rm -f $(PWD)/src/*.d
	rm -f $(PWD)/src/*.pcm
	rm -f $(PWD)/bin/*.exe
	rmdir bin
