

hmmcons_fast.so: hmmcons_fast.c
	g++ -o $@ -O3 -shared -fPIC $^
