import numpy as np
import numpy.linalg as npl

import scipy.integrate as si

import line_profiler

# Calculates the dot product of rows of two matrices
def dot_rows(a, b):
    return np.einsum('ij,ij->i', a, np.conj(b))

# Normalizes an array of vectors (of dimension (N,3)). Note that we are using einsum instead of norm as it's almost twice as fast
def normalize(a):
    return a / np.sqrt(np.einsum('ij,ij->i', a, np.conj(a))).reshape(-1,1)

class OpticalSystem(object):
    def __init__(self, c, Rp, nr):
        # Particle properties
        self.set_particle_index(nr)
        self._c = c
        self.set_particle_radius(Rp)
        
        # Dynamic variables used during raytracing
        self._o = np.array([0, 0, 0])
        self._l = np.array([0, 0, 0])
        
    def set_particle_radius(self, Rp):
        # Make sure that the sphere radius is not zero or negative
        if Rp <= 0:
            raise ValueError("The sphere radius is not valid: {0}".format(Rp))
        self._Rp = Rp
        
    # Set relative index
    def set_particle_index(self, nr):
        if nr > 0:
            self._nr = nr
        else:
            raise ValueError("The sphere refractive index is not valid: {0}".format(nr))
    
    # Calculates the refraction angle given the incidence angle and the relative index of refraction
    # theta can be a 1D numpy array. Then, the return value will be also be a numpy array (Snell's law applied to each element)
    def _snell(self, theta):
        # Raise exception if the angle is out of the [0,pi/2] range
        if np.any((theta < 0) | (theta > np.pi/2)):
            raise ValueError("Incidence angle out of range")
        
        return np.arcsin(1/self._nr * np.sin(theta))
    
    # Calculates the transmission and reflection for a ray with a given incidence angle (th), refraction angle (r), and polarization angle (p) when the relative index of refraction is specified (nr)
    # Also note that the polarization is specified as the normalized power of p-polarization (Pp). Then, the power of the s-polarization is simply (1-Pp).
    # The arguments can be numpy arrays, but then they will have to be 1D and have the same length.
    def _fresnel(self, th, r, Pp):
        nr = self._nr
        
        # Calculate the reflectivities:
        costh = np.cos(th)
        cosr = np.cos(r)
        
        Rs = ((costh - nr*cosr)/(costh + nr*cosr))**2
        Rp = ((cosr - nr*costh)/(cosr + nr*costh))**2
        
        # Calculate the final reflectivity and transmittivity
        R = Rs*(1-Pp) + Rp*Pp
        T = 1 - R
        
        # And return the value
        return (T, R)
    
    def _intersection_angle(self):
        # Make l (director of the line) unitary
        self._l = normalize(self._l)
        ln = self._l
        
        # Note: self._o and self_c should be (N x 3) matrices with N the number of rays considered
        oc = self._o - self._c
        
        # Calculate the discriminant (to see whether there are any solutions)
        
        # The dot products between ln's and oc's because it will be used a lot later
        ln_dot_oc = dot_rows(ln, oc)
        
        # Norms squared of oc's (because this operation is more efficient)
        oc_dot_oc = dot_rows(oc, oc)
        
        D = ln_dot_oc**2 - oc_dot_oc + self._Rp**2
        
        # Discriminant values below zero indicate no intersection, which we will denote by NaN
        D[D < 0] = np.nan
        
        # For convenience, calculate the square root of the determinants, as this will be used later a couple of times
        sqrtD = np.sqrt(D)
        
        # The distance along the line where an intersection occurs (doesn't matter which since this is a sphere)
        d = -ln_dot_oc + sqrtD
        
        # The points at which the intersections occur is x:
        x = self._o + d.reshape(-1,1)*ln
        
        # The vector that points from the center of the sphere to the intersection is r:
        r = x - self._c
        
        # To find the angle of intersecting ray with the normal to the surface, first find the cosine of that angle (absolute value of it)
        c_angles = np.abs(dot_rows(ln, r)/self._Rp)
        
        # Sometimes due to floating-point errors, the value will be slightly higher than 1. The following corrects it:
        c_angles[(c_angles > 1) & (c_angles < 1+1e-8)] = 1
        
        # And finally, return the angle (absolute value)
        return np.arccos(c_angles)
    
    # This function calculates the normalized force (i.e. actual force multiplied by c/(n_1 P)) of a single ray described by a line whose origin is o and whose direction of propagation is l. The sphere of radius R has its center in c and has refractive index nr.
    # Important note: the polarization p is a Jones' vector specified in the lab's coordinate system (e.g. before entering the lens, so that it only has XY components). This vector can be complex. For example, for circular polarization this vector would be (1,i,0), while for linear polarization it is completely real. Its normalization is not important as it is normalized in the code.
    def _ray_force(self, p):
        # Calculate the incidence angle first. NaN values will be passed because they will be filtered later
        th = self._intersection_angle()
        
        ## First we have to determine the coordinate system for the gradient and scattering forces:
        # The scattering force direction, according to Ashkin, 1992, is along the ray propagation direction. Since we have normalized it before in intersection_angle, we don't have to do it again
        dir_scat = self._l
        
        # The gradient force direction (Ashkin, 1992) is orthogonal to the ray propagation direction and lies in the plane formed by the ray and the center of the sphere. For that, we first make a vector that points from the center of the sphere to one of the points in the line and Gram-Schmidt orthogonalize it to make a vector perpendicular to the scattering
        a = self._o - self._c
        
        dir_grad = a - dot_rows(a, dir_scat).reshape(-1,1)*dir_scat
        
        # If the rays pass through the center of the sphere, then dir_grad will be = 0, but this is not a problem as this ray will not exert any gradient force. Then, we can take any direction as dir_grad without any consequence. Otherwise, we have to normalize dir_grad
        dir_grad = normalize(dir_grad)# / npl.norm(dir_grad, axis=1).reshape(-1,1)
        
        # Now remove the undefined values (division by zero)
        dir_grad[np.isnan(dir_grad)] = 0
            
        # The magnitudes of the forces are specified in Ashkin, 1992. First let's calculate some auxiliary quantities:
        # Refraction angles:
        r = self._snell(th)
        
        # Transmission and reflection coefficients
        # Let's calculate the projection of the polarization vector on the incidence plane and the magnitude of that projection
        # First, normalize calculate the norm squared of the polarization for each ray
        pn = dot_rows(p, p)
        
        # Then, use it to calculate the projection of the polarization on the incidence plane
        Pp = (np.abs(dot_rows(p, dir_grad))**2 + np.abs(dot_rows(p, dir_scat))**2)/pn
        
        # Sometimes, the proportion will be slightly bigger than 1 because of floating-point errors. The following corrects it:
        Pp[(Pp > 1) & (Pp < 1+1e-7)] = 1
        
        # Note: if dir_grad is null (when the ray is normal on the sphere), Pp will take some value between 0 and 1, but it won't matter since at normal incidence, Fresnel doesn't depend on the polarization
        
        T, R = self._fresnel(th, r, Pp)
        
        # And finally, the scattering force magnitude:
        # First, some auxiliary arrays (to not compute them twice):
        th2 = 2*th
        r2 = 2*r
        Rsin2th = R*np.sin(th2)
        Rcos2th = R*np.cos(th2)
        denominator = (1 + R**2 + 2*R*np.cos(r2))
        th2_r2 = th2-r2
        Tsq = T**2

        Fs = 1 + Rcos2th - (Tsq * (np.cos(th2_r2) + Rcos2th)) / denominator
        Fs = Fs.reshape(-1,1)
        
        # And the gradient force magnitude:
        Fg = Rsin2th - (Tsq * (np.sin(th2_r2) + Rsin2th)) / denominator
        Fg = Fg.reshape(-1,1)
        
        # And calculate the total force:
        # Note that the sign of Fg is due to a sign error (or maybe misunderstanding?) in Ashkin, 1992
        F = Fs*dir_scat - Fg*dir_grad
        
        # All the non-intersection forces will now be null
        F[np.isnan(F)] = 0
        
        return F
  
# An optical system where all the rays are focused into a single spot (most common arrangement)  
class OpticalSystemSimple(OpticalSystem):
    def __init__(self, c, Rp, nr, Rl, f, p):
        super().__init__(c, Rp, nr)
        
        self._p_single = p
        
        self.set_focal_distance(f)
        self.set_lens_radius(Rl)
        
        self._c = np.array([np.array([0, 0, f]) + c])
        
    def set_focal_distance(self, f):
        if f > 0:
            self._f = f
        else:
            raise ValueError("Invalid focal distance: {0}".format(f))
        
    def set_lens_radius(self, Rl):
        if Rl > 0:
            self._Rl = Rl
        else:
            raise ValueError("Invalid lens radius: {0}".format(Rl))
        
    # Sets the position of the particle relative to the focal spot
    def set_particle_center(self, c):
        self._c = np.array([np.array([0, 0, self._f]) + c])
        
    # Generates the ray directions and origins for a list of r's and th's on the lens (assuming that all the rays are focused in the focal spot)
    def _gen_ray_directions(self, r, th):
        n_rays = len(r)
        
        self._o = np.array([r*np.cos(th), r*np.sin(th), np.zeros(n_rays)]).transpose()
        self._l = np.tile(np.array([0, 0, self._f]), (n_rays, 1)) - self._o
        
        return
    
    # Calculate forces for each and every ray. To be inherited and implemented in children classes
    def _total_ray_force(self, rs, ths):
        pass
    
    # Integrates all the rays, dividing the lens radius by rsteps and the polar angle (2pi) into thsteps
    def integrate(self, rsteps, thsteps):
        rrange = np.linspace(0, self._Rl, rsteps)
        thrange = np.linspace(0, 2*np.pi, thsteps)
        
        # Create the values on which the function will be evaluated
        rs,ths = np.meshgrid(rrange, thrange)
        rs = rs.flatten()
        ths = ths.flatten()
        
        dr = self._Rl/(rsteps-1)
        dth = 2*np.pi/(thsteps-1)
        
        forces = self._total_ray_force(rs, ths)
        Ft = dr*dth*np.sum(forces, axis=0)
        
        return Ft

# A simple system where the intensity on the lens is constant and all the rays are focussed into a single spot
class OpticalSystemSimpleUniform(OpticalSystemSimple):
    def __init__(self, c, Rp, nr, Rl, f, p):
        super().__init__(c, Rp, nr, Rl, f, p)
                
    # Returns the total force by single rays (multiplied by r for polar integration)
    def _total_ray_force(self, r, th):
        n_rays = len(r)
        super()._gen_ray_directions(r, th)
        
        # Make the polarization vectors have the correct dimension
        self._p = np.tile(self._p_single, (n_rays, 1))
        
        F = self._ray_force(self._p)
    
        # The constant factor is to have unit power
        return (r/(np.pi * self._Rl**2)).reshape(-1,1)*F