#include "main.h"

#define dblarg(x) atof(argv[x])

double intens(double r, double th)
{
    return exp(-pow(r, 2));
}

int main(int argc, char* argv[])
{
    // Here come the arguments for the simulation
    
    // The lens radius
    const double lens_r = dblarg(1);
    
    // The lens focal distance
    const double df = dblarg(2);
    
    // Radius, external index and internal index of the sphere
    const double sph_r = dblarg(3);
    const double sph_ne = dblarg(4);
    const double sph_n = dblarg(5);
    
    // Quantity of steps to divide the r and theta (integration variables) into.
    const double r_steps = dblarg(6);
    const double th_steps = dblarg(7);
    
    // The limits of integration for x (radial position of particle) IN TERMS OF
    // PARTICLE RADIUS. The particle position is relative to the focus.
    const double x_final = dblarg(8)*sph_r;
    const double x_steps = dblarg(9);
    
    // The limits of integration for z (axial position of particle) IN TERMS OF
    // PARTICLE RADIUS. The particle position is relative to the focus.
    const double z_init = dblarg(10)*sph_r;
    const double z_final = dblarg(11)*sph_r;
    const double z_steps = dblarg(12);
    
    // Calculate the differentials
    const double dr = 1/r_steps;
    const double dth = 2*M_PI/th_steps;
    
    const double dx = x_final/x_steps;
    const double dz = z_final/z_steps;
    
    Sphere s = Sphere();
    Lens l = Lens();
    
    l.set_df(df);
    l.set_intensity(&intens);
    l.set_radius(lens_r);
    
    s.set_r(sph_r);
    s.set_n(sph_n);
    s.set_ne(sph_ne);
    
    // Begin the simulation
    
    Vector3d force;
    
    for (double z=0; z <= z_final; z += dz)
    {
        for (double x=0; x <= x_final; x += dx)
        {
            force = Vector3d(0,0,0);
            l.set_lens_pos(Vector3d(-x,0,df-z));
            for (double r=0; r<=1; r+= dr)
            {
                for (double th = 0; th <= 2*M_PI; th+= dth)
                {
                    force += s.get_force(l.get_ray(r, th))*dr*dth;
                }
            }
            printf("%e %e %e %e\n", x, z, force[0], force[2]);
        }
    }
    
    return 0;
}
