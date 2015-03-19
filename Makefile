LIBS=-lrt ./htslib/libhts.a -lz
CFLAGS=-fopenmp -O3

SRC=nanopolish.cpp \
    nanopolish_consensus.cpp \
    nanopolish_khmm_parameters.cpp \
    nanopolish_klcs.cpp \
    nanopolish_common.cpp \
    nanopolish_khmm.cpp \
    nanopolish_profile_hmm.cpp \
    nanopolish_anchor.cpp

nanopolish: $(SRC)
	g++ -o $@ $(CFLAGS) -fPIC $^ $(LIBS)

clean:
	rm nanopolish
