#
LIBS=-lrt ./htslib/libhts.a -lz
CPPFLAGS=-fopenmp -O3
PROGRAM=nanopolish

default: $(PROGRAM)

# Source files
SRC=nanopolish.cpp \
    nanopolish_consensus.cpp \
    nanopolish_khmm_parameters.cpp \
    nanopolish_klcs.cpp \
    nanopolish_common.cpp \
    nanopolish_khmm.cpp \
    nanopolish_profile_hmm.cpp \
    nanopolish_anchor.cpp

# Automatically generated object names
OBJ=$(SRC:.cpp=.o)

# Generate dependencies
depend: .depend

.depend: $(SRC)
	rm -f ./.depend
	g++ $(CPPFLAGS) -MM $^ > ./.depend;

include .depend

# Compile objects
.cpp.o:
	g++ -c $(CPPFLAGS) -fPIC $<

# Link executable
$(PROGRAM): $(OBJ)
	g++ -o $@ $(CPPFLAGS) -fPIC $^ $(LIBS)

clean:
	rm nanopolish
