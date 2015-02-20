LIBS=-lrt
CFLAGS=-fopenmp

libnanopolish.so: nanopolish.cpp nanopolish_khmm_parameters.cpp
	g++ -o $@ $(CFLAGS) -O3 -shared -fPIC $^ $(LIBS)
