

hmmcons_fast.so: hmmcons_fast.c
	g++ -o $@ -g -O2 -shared -fPIC $^ -lrt
