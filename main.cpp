#include "main.h"

#define dblarg(x) atof(argv[x])

// Since we accept normalized r, the beam width doesn't appear.
// The total laser power is unity to independize the calculation of force 
double intens(double r, double th)
{
	const double I0 = 2/M_PI;
    return I0*exp(-2*pow(r, 2));
}

int main(int argc, char* argv[])
{
    // Here come the arguments for the simulation
    
    // The lens radius will be unity
    // The lens focal distance is calculated from NA
    const double NA = dblarg(1);
    
    // The refractive index of the medium
    const double ne = dblarg(2);
    
    // Now we calculate the focal distance
    const double df = sqrt(pow(ne, 2) - pow(NA, 2))/NA;
    
    // Index of refraction of the sphere (relative)
    const double sph_n = dblarg(3);
    
    // Quantity of steps to divide the r and theta (integration variables) into.
    const double r_steps = dblarg(4);
    const double th_steps = dblarg(5);
    
    // The limits of integration for x (radial position of particle) IN TERMS OF
    // PARTICLE RADIUS. The particle position is relative to the focus.
    const double x_init = dblarg(6);
    const double x_final = dblarg(7);
    const double x_steps = dblarg(8);
    
    // The limits of integration for z (axial position of particle) IN TERMS OF
    // PARTICLE RADIUS. The particle position is relative to the focus.
    const double y_init = dblarg(9);
    const double y_final = dblarg(10);
    const double y_steps = dblarg(11);
    
    // The limits of integration for z (axial position of particle) IN TERMS OF
    // PARTICLE RADIUS. The particle position is relative to the focus.
    const double z_init = dblarg(12);
    const double z_final = dblarg(13);
    const double z_steps = dblarg(14);
    
    // Calculate the differentials
    const double dr = 1/r_steps;
    const double dth = 2*M_PI/th_steps;
    
    const double dx = (x_final-x_init)/x_steps;
    const double dy = (y_final-y_init)/y_steps;
    const double dz = (z_final-z_init)/z_steps;
    
    Sphere s = Sphere();
    Lens l = Lens();
    
    l.set_df(df);
    l.set_intensity(&intens);
    l.set_radius(1);
    
    s.set_r(1);
    s.set_n(sph_n*ne);
    s.set_ne(ne);
    
    // Begin the simulation
    
    Vector3d force;
    
    // Add tolerance to the comparations so that the rounding errors don't
    // influence the quantity of cycles.
    
    double tolx = (x_final - x_init)/1e4;
    double toly = (y_final - y_init)/1e4;
    double tolz = (z_final - z_init)/1e4;
    
    for (double x=x_init; x <= x_final+tolx; x += dx)
    {
        for (double y=y_init; y <= y_final+toly; y += dy)
        {
        	for (double z=z_init; z <= z_final+tolz; z += dz)
     	   	{
		        force = Vector3d(0,0,0);
		        l.set_lens_pos(Vector3d(-x,0,df-z));
		        
		        // Get the force and independize itforce from the external 
		        // refractive index and the laser power
		        force = s.get_total_force(l, dr, dth) * c / ne;
		        
		        printf("%e %e %e %e %e %e\n", x, y, z, force[0], force[1], force[2]);
		        if (z_init == z_final) break;
	        }
	        printf("\n");
	        if (y_init == y_final) break;
        }
        printf("\n");
        if (x_init == x_final) break;
    }
    
    return 0;
}
