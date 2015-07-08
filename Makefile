#

# Sub directories containing source code, except for the main programs
SUBDIRS := src src/hmm

#
# Set up libraries and paths
#

# Default to automatically installing hdf5
H5_LIB=./lib/libhdf5.a -ldl
H5_INCLUDE=-I./include
HTS_LIB=./htslib/libhts.a


HTS_INCLUDE=-I./htslib
FAST5_INCLUDE=-I./fast5
NP_INCLUDE=$(addprefix -I./, $(SUBDIRS))

LIBS=-lrt $(HTS_LIB) -lz $(H5_LIB)
CPPFLAGS=-fopenmp -O3 -std=c++11 -g $(H5_INCLUDE) $(HTS_INCLUDE) $(FAST5_INCLUDE) $(NP_INCLUDE)
CFLAGS=-O3

PROGRAM=nanopolish
TEST_PROGRAM=nanopolish_test

CXX=g++
CC=gcc

all: $(PROGRAM) $(TEST_PROGRAM)

#
# Build libhts
#
htslib/libhts.a:
	cd htslib; make

#
# Automatically install HDF5 dependency if requested by user
#
lib/libhdf5.a:
	wget https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz
	tar -xzf hdf5-1.8.14.tar.gz
	cd hdf5-1.8.14; ./configure --enable-threadsafe --prefix=`pwd`/..; make; make install

# Overwrite H5 variables to use system-wide version
.PHONY: libhdf5.system
libhdf5.system:
	$(eval H5_LIB=-lhdf5)
	$(eval H5_INCLUDE=)

# Find the source files by searching subdirectories
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard $(dir)/*.c))

EXE_SRC=src/main/nanopolish.cpp src/test/nanopolish_test.cpp

# Automatically generated object names
CPP_OBJ=$(CPP_SRC:.cpp=.o)
C_OBJ=$(C_SRC:.c=.o)

# Generate dependencies
PHONY=depend
depend: .depend

.depend: $(CPP_SRC) $(C_SRC) $(EXE_SRC) $(H5_LIB)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS) -MM $(CPP_SRC) $(C_SRC) > ./.depend;

include .depend

# Compile objects
.cpp.o:
	$(CXX) -o $@ -c $(CPPFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) -fPIC $<


# Link main executable
$(PROGRAM): src/main/nanopolish.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB)
	$(CXX) -o $@ $(CPPFLAGS) -fPIC $^ $(LIBS)

# Link test executable
$(TEST_PROGRAM): src/test/nanopolish_test.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB)
	$(CXX) -o $@ $(CPPFLAGS) -fPIC $^ $(LIBS)

test: $(TEST_PROGRAM)
	./$(TEST_PROGRAM)

clean:
	rm nanopolish nanopolish_test $(CPP_OBJ) $(C_OBJ) src/main/nanopolish.o src/test/nanopolish_test.o
