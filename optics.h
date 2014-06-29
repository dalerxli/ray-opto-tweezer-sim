#include <iostream>
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;
using namespace std;

//#define DISABLE_REFLECTION

// Auxiliary functions

double snell(double thi, double nr);
double R_Fresnel(double thi, double nr);

// Constants
const double c = 299792458.0; // Light speed

// Represents a light ray with a specified power.
class Ray
{
    private:
        // The first point through which the ray passes
        Vector3d origin;
        
        // A unitary vector pointing in the direction of propagation
        Vector3d direction;
        
        // The optical power of the ray
        double power;
    
    public:
        // Set the direction of propagation and the origin. The direction is 
        // normalized
        void set_vect(Vector3d orig, Vector3d p2);
        
        // Set the power
        void set_power(double pow){this->power = pow;}
        
        // Get the reference to the origin
        Vector3d& get_orig(void) {return this->origin;}
        
        // Get the reference to the direction of propagation (normalized)
        Vector3d& get_dir(void) {return this->direction;}
        
        // Get the power
        double& get_pow(void) {return this->power;}
};

// The focusing lens (e.g. microscope objective). Considered to be thin. The 
// focal point is considered to be below the lens at (x0, y0, z0-df), where 
// (x0,y0,z0) is the position of the center of the lens and df is the focal 
// distance.
class Lens
{
    private:
        double radius;
        double df;
        Vector3d lens_pos;
        Vector3d focus_pos;
        
        // A ray that is currently "active". This avoids the creation of many
        // Ray instances.
        Ray current_ray;
        
        // A function f(r, th) depending on _normalized_ radial distance and 
        // polar angle that gives the intensity in that point of the lens
        double (*intensity_f)(double r, double th);
        
    public:
        void set_lens_pos(Vector3d pos);
        void set_df(double df);
        double get_df(){return this->df;}
        void set_radius(double r){this->radius = r;}
        double get_radius(){return this->radius;}
        void set_intensity(double (*in_f)(double, double))
        {
            this->intensity_f = in_f;
        }
        
        // Get a ray from a parametrized point in the sphere. The radius is 
        // normalized. Accepts the position in form (r, theta) and returns the ray
        // with the intensity corrected for polar integration.
        Ray& get_ray(double r, double th);
};

// Represents a dielectric sphere with refractive index n, external refractive 
// index ne, radius r. The position is always (0,0,0)

class Sphere
{
    private:
        double n;
        double ne;
        double nr; // Relative index
        double r;
        
        // A force that is currently "active". This avoids the creation of many
        // Vector3d instances.
        Vector3d current_force;
        
    public:
        void set_n(double n)
        {
            this->n = n;
            this->nr = n/this->ne;
        }
        
        void set_ne(double ne)
        {
            this->ne = ne;
            this->nr = this->n/ne;
        }
        
        void set_r(double r) {this->r = r;}
        
        // Force due to one ray
        Vector3d& get_force(Ray& ray);

		// Force due to all the rays with the specified steps
		Vector3d get_total_force(Lens& l, double dr, double dth);
};
