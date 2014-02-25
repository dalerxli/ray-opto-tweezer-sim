all:
	g++ -I /usr/include/eigen3/ -g main.cpp optics.cpp -O2 -march=atom -mtune=atom -o app
