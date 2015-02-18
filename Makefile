
#LIBS=-lrt
#CFLAGS=-fopenmp

hmmcons_fast.so: hmmcons_fast.cpp hmmcons_khmm_parameters.cpp
	g++ -o $@ $(CFLAGS) -O3 -shared -fPIC $^ $(LIBS)
