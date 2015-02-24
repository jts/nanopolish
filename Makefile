LIBS=-lrt
CFLAGS=-fopenmp -O3
SRC=nanopolish.cpp nanopolish_khmm_parameters.cpp nanopolish_klcs.cpp nanopolish_common.cpp nanopolish_khmm.cpp nanopolish_profile_hmm.cpp

libnanopolish.so: $(SRC)
	g++ -o $@ $(CFLAGS) -shared -fPIC $^ $(LIBS)

clean:
	rm libnanopolish.so
