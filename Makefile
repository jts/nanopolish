#

#
# Set up libraries and paths
#

# Default to using the system-wide libhdf
H5_LIB=-lhdf5
H5_INCLUDE=

LIBS=-lrt ./htslib/libhts.a -lz $(H5_LIB)
CPPFLAGS=-fopenmp -O3 -std=c++11 -g $(H5_INCLUDE)
PROGRAM=nanopolish
CXX=g++

all: $(PROGRAM)

#
# Build libhts
#
htslib/libhts.a:
	cd htslib; make

#
# Automatically install HDF5 dependency if requested by user
#
lib/libhdf5.a:
	wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.14.tar.gz
	tar -xzf hdf5-1.8.14.tar.gz
	cd hdf5-1.8.14; ./configure --prefix=`pwd`/..; make; make install

# Overwrite H5 variables to put to local version
.PHONY: libhdf5.install
libhdf5.install: lib/libhdf5.a
	$(eval H5_LIB=./lib/libhdf5.a -ldl)
	$(eval H5_INCLUDE=-I./include)

# Source files
SRC=nanopolish.cpp \
    nanopolish_consensus.cpp \
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
    logsum.cpp

# Automatically generated object names
OBJ=$(SRC:.cpp=.o)

# Generate dependencies
depend: .depend

.depend: $(SRC)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS) -MM $^ > ./.depend;

include .depend

# Compile objects
.cpp.o:
	$(CXX) -c $(CPPFLAGS) -fPIC $<

# Link executable
$(PROGRAM): $(OBJ) htslib/libhts.a
	$(CXX) -o $@ $(CPPFLAGS) -fPIC $^ $(LIBS)

clean:
	rm nanopolish *.o
