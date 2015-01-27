

hmmcons_fast.so: hmmcons_fast.c
	g++ -Wall -o $@ -g -O2 -shared -fPIC $^ -lrt
