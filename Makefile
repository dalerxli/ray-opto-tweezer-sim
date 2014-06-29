all:
	g++ -I /usr/include/eigen3/ -g main.cpp optics.cpp -march=native -O2 -o app
test:
	g++ -I /usr/include/eigen3/ -g test_ray_q.cpp optics.cpp -march=native -O2 -o test_q
