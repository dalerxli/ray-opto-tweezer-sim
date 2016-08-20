all:
	g++ -I /usr/include/eigen3/ -g main.cpp optics.cpp cuba/libcuba.so -march=native -O2 -o app
