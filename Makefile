all:
	g++ -I /usr/include/eigen3/ -g main.cpp optics.cpp -O2 -march=atom -mtune=atom -o app
test:
	g++ -I /usr/include/eigen3/ -g test_ray_q.cpp optics.cpp -O2 -march=atom -mtune=atom -o test_q
