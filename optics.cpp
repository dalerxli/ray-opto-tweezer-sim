#include "optics.h"

// Auxiliary functions

// Snell law
double snell(double thi, double nr)
{
    return asin(sin(thi)/nr);
}

// Fresnel reflection coefficient for circularly polarized light. nr is the
// relative index

// TODO: add polarization.
double R_Fresnel(double thi, double nr)
{
    #ifdef DISABLE_REFLECTION
    return 0;
    #endif
    
    double thr = snell(thi, nr);
    double R = pow((cos(thi) - nr*cos(thr))/(cos(thi) + nr*cos(thr)), 2);
    R += pow((cos(thr) - nr*cos(thi))/(cos(thr) + nr*cos(thi)), 2);
    
    return R/2;
}

// Ray methods

void Ray::set_vect(Vector3d orig, Vector3d p2)
{
    this->origin = orig;
    this->direction = (p2-orig).normalized();
    
    return;
}

// Lens methods
void Lens::set_lens_pos(Vector3d pos)
{
    this->lens_pos = pos;
    this->focus_pos = pos - Vector3d(0,0,this->df);
    return;
}

void Lens::set_df(double df)
{
    this->df = df;
    this->focus_pos = this->lens_pos - Vector3d(0,0,df);
    return;
}

Ray& Lens::get_ray(double rn, double th)
{
    double r = rn * this->radius;
    double P = rn * this->intensity_f(rn, th);
    Vector3d polar = r*Vector3d(cos(th), sin(th), 0);
    
    this->current_ray.set_vect(this->lens_pos+polar, this->focus_pos);
    this->current_ray.set_power(r*P);
    
    return this->current_ray;
}

// Sphere methods
Vector3d& Sphere::get_force(Ray& ray)
{
    Vector3d l = ray.get_dir();
    Vector3d x0 = ray.get_orig();
    
    double discr = pow(l.dot(x0),2) - x0.dot(x0) + pow(this->r, 2);
    
    // In the following case there is either one tangential intersection or
    // no intersection at all.
    if (discr <= 0)
    {
        this->current_force = Vector3d(0,0,0);
        return this->current_force;
    }
        
    // Cosine of the incidence angle
    double costh = sqrt(discr)/this->r;
    
    Vector3d lperp = Vector3d(0,0,0);
    
    if (costh < 1.0)
    {
        // The normal vector (equal to the position of intersection because
        // the sphere is centered). We take the smallest value because we need
        // the first intersection. Note that the dot product is negative
        double d = -l.dot(x0) - sqrt(discr);
        Vector3d normal = (x0 + d*l).normalized();
        
        // Use the Gram-Schmidt ortogonalization to find the orthogonal vector
        lperp = (normal - normal.dot(l)*l).normalized();
    }
    
    // Sometimes the number is marginally greater than 1 (due to precision 
    // limits)
    else if (costh > 1.0) costh = 1;
    
    // Now we calculate the force
    double thi = acos(costh);
    double thr = snell(thi, this->nr);
    
    // Force coefficient
    double k = ray.get_pow() * this->ne / c;
    
    // Fresnel coefficients
    double R = R_Fresnel(thi, this->nr);
    double T = 1 - R;
    
    // Double angles
    double thi2 = thi * 2;
    double thr2 = thr * 2;
    
    // An auxiliary coefficient
    double al = 1 + pow(R, 2) + 2*R*cos(thr2);
    
    // Now the components of the force
    double Fs = k*(1+R*cos(thi2) - pow(T, 2)*(cos(thi2 - thr2) + R*cos(thi2))/al);
    double Fg = -k*(R*sin(thi2) - pow(T, 2)*(sin(thi2 - thr2) + R*sin(thi2))/al);
    
    this->current_force = Fs*l + Fg*lperp;
    
    if (isnan(this->current_force[0]))
    {
        this->current_force = Vector3d(0,0,0);
    }
    return this->current_force;
}
