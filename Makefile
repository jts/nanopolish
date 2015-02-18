
#LIBS=-lrt
#CFLAGS=-fopenmp

hmmcons_fast.so: hmmcons_fast.c hmmcons_khmm_parameters.c
	g++ -o $@ $(CFLAGS) -O3 -shared -fPIC $^ $(LIBS)
