

hmmcons_fast.so: hmmcons_fast.c hmmcons_khmm_parameters.c
	g++ -o $@ -g -fopenmp -O2 -shared -fPIC $^ -lrt
