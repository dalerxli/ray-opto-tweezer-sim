#include "main.h"
#define dblarg(x) atof(argv[x])

// This file tests the Q quantities for rays with different angles of incidence
// They have to coincide more or less with those of Ashkin

// Since we accept normalized r, the beam width doesn't appear 
double intens(double r, double th)
{
    return exp(-2*pow(r, 2));
}

int main(int argc, char* argv[])
{
    Sphere s = Sphere();
    
    // Q doesn't depend on particle's radius
    s.set_r(1);
    s.set_n(dblarg(1));
    
    const double ne = dblarg(2);
    s.set_ne(ne);
    
    Ray ray = Ray();
    ray.set_power(1);
    
    Vector3d p2 = Vector3d(0,0,1);
    
    const double dth = (M_PI/2)/(90*2);
    for (double th=0; th <= M_PI/2; th += dth)
    {
    	Vector3d orig = Vector3d(0, sin(th), cos(th));
    	ray.set_vect(orig, p2);
    	
    	Vector3d f = s.get_force(ray);
    	double Q = f.norm()*c/(ne*ray.get_pow());
    	
    	printf("%e %e\n", th/M_PI*180, Q);
    }
    
    return 0;
}
