#
LIBS=-lrt ./htslib/libhts.a -lz -lhdf5
CPPFLAGS=-fopenmp -O3 -std=c++11 -g
PROGRAM=nanopolish
CXX=g++

default: $(PROGRAM)

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
$(PROGRAM): $(OBJ)
	$(CXX) -o $@ $(CPPFLAGS) -fPIC $^ $(LIBS)

clean:
	rm nanopolish *.o
