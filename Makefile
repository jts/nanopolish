#

#
# Set up libraries and paths
#

# Default to automatically installing hdf5
H5_LIB=./lib/libhdf5.a -ldl
H5_INCLUDE=-I./include
HTS_LIB=./htslib/libhts.a

LIBS=-lrt $(HTS_LIB) -lz $(H5_LIB)
CPPFLAGS=-fopenmp -O3 -std=c++11 -g $(H5_INCLUDE)
PROGRAM=nanopolish
TEST_PROGRAM=nanopolish_test

CXX=g++

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

# Source files, except for the main programs
SRC=nanopolish_consensus.cpp \
    nanopolish_khmm_parameters.cpp \
    nanopolish_klcs.cpp \
    nanopolish_common.cpp \
    nanopolish_profile_hmm.cpp \
    nanopolish_anchor.cpp \
    nanopolish_fast5_map.cpp \
    nanopolish_poremodel.cpp \
    nanopolish_squiggle_read.cpp \
    nanopolish_eventalign.cpp \
    nanopolish_getmodel.cpp \
    nanopolish_iupac.cpp \
    logsum.cpp

EXE_SRC=nanopolish.cpp nanopolish_test.cpp

# Automatically generated object names
OBJ=$(SRC:.cpp=.o)

# Generate dependencies
PHONY=depend
depend: .depend

.depend: $(SRC) $(EXE_SRC) $(H5_LIB)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS) -MM $(SRC) > ./.depend;

include .depend

# Compile objects
.cpp.o:
	$(CXX) -c $(CPPFLAGS) -fPIC $<

# Link main executable
$(PROGRAM): nanopolish.o $(OBJ) $(HTS_LIB) $(H5_LIB)
	$(CXX) -o $@ $(CPPFLAGS) -fPIC $^ $(LIBS)

# Link test executable
$(TEST_PROGRAM): nanopolish_test.o $(OBJ) $(HTS_LIB) $(H5_LIB)
	$(CXX) -o $@ $(CPPFLAGS) -fPIC $^ $(LIBS)

test: $(TEST_PROGRAM)
	./$(TEST_PROGRAM)

clean:
	rm nanopolish *.o
