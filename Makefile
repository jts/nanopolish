#

# Sub directories containing source code, except for the main programs
SUBDIRS := src src/hmm src/thirdparty src/thirdparty/scrappie src/common src/alignment src/pore_model

#
# Set libraries, paths, flags and options
#

#Basic flags every build needs
LIBS=-lz
CXXFLAGS ?= -g -O3
CXXFLAGS += -std=c++11 -fopenmp -fsigned-char
CFLAGS ?= -O3 -std=c99
CXX ?= g++
CC ?= gcc

# Change the value of HDF5, EIGEN, or HTS below to any value to disable compilation of bundled code
HDF5?=install
EIGEN?=install
HTS?=install

# Check operating system, OSX doesn't have -lrt
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LIBS += -lrt
endif

# Default to automatically installing hdf5
ifeq ($(HDF5), install)
    H5_LIB=./lib/libhdf5.a
    H5_INCLUDE=-I./include
    LIBS += -ldl
else
    # Use system-wide hdf5
    H5_LIB=
    H5_INCLUDE=
    LIBS += -lhdf5
endif

# Default to automatically installing EIGEN
ifeq ($(EIGEN), install)
    EIGEN_CHECK=eigen/INSTALL
else
    # Use system-wide eigen
    EIGEN_CHECK=
endif

# Default to build and link the libhts submodule
ifeq ($(HTS), install)
    HTS_LIB=./htslib/libhts.a
    HTS_INCLUDE=-I./htslib
else
    # Use system-wide htslib
    HTS_LIB=
    HTS_INCLUDE=
    LIBS += -lhts
endif

# Include the header-only fast5 library
FAST5_INCLUDE=-I./fast5/include

# Include the header-only eigen library
EIGEN_INCLUDE=-I./eigen/

# Include the src subdirectories
NP_INCLUDE=$(addprefix -I./, $(SUBDIRS))

# Add include flags
CPPFLAGS += $(H5_INCLUDE) $(HTS_INCLUDE) $(FAST5_INCLUDE) $(NP_INCLUDE) $(EIGEN_INCLUDE)

# Main programs to build
PROGRAM=nanopolish
TEST_PROGRAM=nanopolish_test

all: $(PROGRAM) $(TEST_PROGRAM)

#
# Build libhts
#
htslib/libhts.a:
	cd htslib && make || exit 255

#
# If this library is a dependency the user wants HDF5 to be downloaded and built.
#
lib/libhdf5.a:
	if [ ! -e hdf5-1.8.14.tar.gz ]; then wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz; fi
	tar -xzf hdf5-1.8.14.tar.gz || exit 255
	cd hdf5-1.8.14 && ./configure --enable-threadsafe --prefix=`pwd`/.. && make && make install


# Download and install eigen if not already downloaded
eigen/INSTALL:
	if [ ! -e 3.2.5.tar.bz2 ]; then wget http://bitbucket.org/eigen/eigen/get/3.2.5.tar.bz2; fi
	tar -xjf 3.2.5.tar.bz2 || exit 255
	mv eigen-eigen-bdd17ee3b1b3 eigen || exit 255

#
# Source files
#

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

.depend: $(CPP_SRC) $(C_SRC) $(EXE_SRC) $(H5_LIB) $(EIGEN_CHECK)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $(CPP_SRC) $(C_SRC) > ./.depend;

include .depend

# Compile objects
.cpp.o:
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) -fPIC $<

.c.o:
	$(CC) -o $@ -c $(CFLAGS) $(CPPFLAGS) $(H5_INCLUDE) -fPIC $<

# Link main executable
$(PROGRAM): src/main/nanopolish.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB) $(EIGEN_CHECK)
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $< $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB) $(LIBS) $(LDFLAGS)

# Link test executable
$(TEST_PROGRAM): src/test/nanopolish_test.o $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB)
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $< $(CPP_OBJ) $(C_OBJ) $(HTS_LIB) $(H5_LIB) $(LIBS) $(LDFLAGS)

test: $(TEST_PROGRAM)
	./$(TEST_PROGRAM)

clean:
	rm -f $(PROGRAM) $(TEST_PROGRAM) $(CPP_OBJ) $(C_OBJ) src/main/nanopolish.o src/test/nanopolish_test.o
