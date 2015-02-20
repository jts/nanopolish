LIBS=-lrt
CFLAGS=-fopenmp

libnanopolish.so: nanopolish.cpp nanopolish_khmm_parameters.cpp nanopolish_klcs.cpp
	g++ -o $@ $(CFLAGS) -O3 -shared -fPIC $^ $(LIBS)
