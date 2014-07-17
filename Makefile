all:
	g++ -I /usr/include/eigen3/ -g main.cpp optics.cpp -march=native -O2 -o app
