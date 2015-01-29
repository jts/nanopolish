

hmmcons_fast.so: hmmcons_fast.c hmmcons_khmm_parameters.c
	g++ -o $@ -g -O2 -shared -fPIC $^ -lrt
