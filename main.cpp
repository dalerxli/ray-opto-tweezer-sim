#include "main.h"

#define dblarg(x) atof(argv[x])

// Since we accept normalized r, the beam width doesn't appear.
// The total laser power is unity to independize the calculation of force 
double intens(double r, double th)
{
	const double I0 = 8/M_PI;
    return I0*exp(-2*pow(2*r, 2));
}

int main(int argc, char* argv[])
{
    // Here come the arguments for the simulation
    
    // The focal distance of the lens
    const double df = dblarg(1);
    
    // The radius of the lens (i.e. the radius of the incoming beam)
    const double l_r = 2*dblarg(2);
    
    // The refractive index of the medium
    const double ne = dblarg(3);
    
    // Index of refraction of the sphere (relative)
    const double sph_n = dblarg(4);
    
    // Quantity of steps to divide the r and theta (integration variables) into.
    const double r_steps = dblarg(5);
    const double th_steps = dblarg(6);
    
    // The limits of integration for x (radial position of particle)
    // The particle position is relative to the focus.
    const double x_init = dblarg(7);
    const double x_final = dblarg(8);
    const double x_steps = dblarg(9);
    
    // The limits of integration for z (axial position of particle)
    // The particle position is relative to the focus.
    const double y_init = dblarg(10);
    const double y_final = dblarg(11);
    const double y_steps = dblarg(12);
    
    // The limits of integration for z (axial position of particle)
    // The particle position is relative to the focus.
    const double z_init = dblarg(13);
    const double z_final = dblarg(14);
    const double z_steps = dblarg(15);
    
    // Calculate forces for a double trap? 1 for true, 0 for false
    unsigned int double_trap = (int) dblarg(16);
    if (double_trap != 0) double_trap = 1;
    
    // Additional parameters for the Gaussian beam
    
    // The sphere radius
    const double s_r = dblarg(17);
    
    // The wavelength of the light
    const double lam = dblarg(18);
    
    // Calculate the differentials
    const double dr = 1/r_steps;
    const double dth = 2*M_PI/th_steps;
    
    const double dx = (x_final-x_init)/x_steps;
    const double dy = (y_final-y_init)/y_steps;
    const double dz = (z_final-z_init)/z_steps;
    
    Sphere s = Sphere();
    Lens l = Lens();
    
    //l.set_df(df);
    l.set_intensity(&intens);
    l.set_radius(l_r);
    
    s.set_r(s_r);
    s.set_n(sph_n*ne);
    s.set_ne(ne);
    
    // Begin the simulation
    
    Vector3d force;
    
    // Add tolerance to the comparations so that the rounding errors don't
    // influence the quantity of cycles.
    
    double tolx = (x_final - x_init)/1e4;
    double toly = (y_final - y_init)/1e4;
    double tolz = (z_final - z_init)/1e4;
    
    // Separator after iterating over z
    const char sepz = ((z_init == z_final) ? ' ' : '\n');
    
    // Separator after iterating over y
    const char sepy = ((y_init == y_final) ? ' ' : '\n');
    
    for (double x=x_init; x <= x_final+tolx; x += dx)
    {
        for (double y=y_init; y <= y_final+toly; y += dy)
        {
        	for (double z=z_init; z <= z_final+tolz; z += dz)
     	   	{
     	   	    // Calculate the parameters for the Gaussian beam
     	   	    double th = atan(l_r/df);
     	   	    double w0 = lam/(M_PI*th);
     	   	    double zr = w0/th;
     	   	    
     	   	    // The radius of curvature (with the fix for the infinite
     	   	    // radius of curvature at beam waist)
     	   	    double R = 1e7;
     	   	    if (fabs(z) >= 1e-12)
     	   	        R = fabs(z*(1 + pow(zr/z, 2)));
     	   	    
     	   	    // The beam width
     	   	    double w = w0*sqrt(1 + pow(z/zr, 2));
     	   	    
     	   	    // Finally, we set the optical parameters corrected to
     	   	    // simulate the Gaussian beam
     	   	    
     	   	    // The effective focal length
     	   	    double df_ef = R*l_r/w;
     	   	    l.set_df(df_ef);
     	   	    
     	   	    // The effective distance to the focal point (to have the
     	   	    // correct beam width in the simulation). The weird expression
     	   	    // is needed to copy the sign of z.
     	   	    double z_ef = ((z > 0) - (z < 0)) *  R;
     	   	    
     	   	    //printf("\n%e, %e\n", z_ef, R);
     	   	    
		        force = Vector3d(0,0,0);
		        
		        // Calculate the forces for a single or a double trap
		        for(int i=0; i <= double_trap; i++)
		        {
		        	Vector3d temp = Vector3d(0,0,0);
		        	
		        	// Note that we are using the corrected values - R and z_ef
				    l.set_lens_pos(Vector3d(-x,0, df_ef - z_ef * pow(-1, i)));
				    
				    // Get the force and independize it from the external 
				    // refractive index and the laser power
				    temp = s.get_total_force(l, dr, dth) * c / ne;
				    temp = temp * pow(0.5, double_trap);
				    temp = temp.cwiseProduct(Vector3d(1,1,pow(-1, i)));
				    
				    force += temp;
		        }
		        
		        printf("%e %e %e %e %e %e\n", x, y, z, force[0], force[1], force[2]);
		        if (z_init == z_final) break;
	        }
	        
	        if (y_init == y_final) break;
	        printf("%c", sepz);
        }
        
        if (x_init == x_final) break;
        printf("%c", sepy);
    }
    
    return 0;
}
