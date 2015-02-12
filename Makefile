

hmmcons_fast.so: hmmcons_fast.c hmmcons_khmm_parameters.c
	g++ -o $@ -fopenmp -O3 -shared -fPIC $^ -lrt
