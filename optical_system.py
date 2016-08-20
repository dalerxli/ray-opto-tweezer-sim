import numpy as np
import numpy.linalg as npl

import scipy.integrate as si

import logging

logging.basicConfig(level=logging.INFO)
lgg = logging.getLogger("intersection_module")

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
        if not 0 <= theta <= np.pi/2:
            raise ValueError("Incidence angle out of range: {0}".format(theta))
        
        return np.arcsin(1/self._nr * np.sin(theta))
    
    # Calculates the transmission and reflection for a ray with a given incidence angle (th), refraction angle (r), and polarization angle (p) when the relative index of refraction is specified (nr)
    # Also note that the polarization is specified as the normalized power of p-polarization (Pp). Then, the power of the s-polarization is simply (1-Pp).
    # The arguments can be numpy arrays, but then they will have to be 1D and have the same length.
    def _fresnel(self, th, r, Pp):
        nr = self._nr
        
        # Calculate the reflectivities:
        # #TODO avoid calling cos(angle) multiple times for the same angles
        Rs = ((np.cos(th) - nr*np.cos(r))/(np.cos(th) + nr*np.cos(r)))**2
        Rp = ((np.cos(r) - nr*np.cos(th))/(np.cos(r) + nr*np.cos(th)))**2
        
        # Calculate the final reflectivity and transmittivity
        R = Rs*(1-Pp) + Rp*Pp
        T = 1 - R
        
        # And return the value
        return (T, R)
    
    def _intersection_angle(self):
        # Make l (director of the line) unitary
        ln = self._l/npl.norm(self._l, axis=1).reshape(-1,1)
        
        # Note: self._o and self_c should be (N x 3) matrices with N the number of rays considered
        oc = self._o - self._c
        
        # Calculate the discriminant (to see whether there are any solutions)
        
        # The dot products between ln's and oc's because it will be used a lot later
        ln_dot_oc = np.einsum('ij,ij->i', ln, oc)
        
        # Norms squared of oc's (because this operation is more efficient)
        oc_dot_oc = np.einsum('ij,ij->i', oc, oc)
        
        D = ln_dot_oc**2 - oc_dot_oc + self._Rp**2
        
        # For convenience, calculate the square root of the determinants, as this will be used later a couple of times
        sqrtD = np.sqrt(D)
        
        lgg.debug(D)
        
        # If there are no solutions, then there is nothing else to do
        # TODO: redo this for numpy arrays
        if D < 0:
            return np.nan
        
        # Otherwise, calculate the distance along the line where an intersection occurs (doesn't matter which since this is a sphere)
        d = -ln_dot_oc + sqrtD
        lgg.debug(d)
        
        # The point at which the intersection occurs is x:
        x = self._o + d*ln
        
        # If x is zero (which would be a problem when calculating its inverse norm), switch to the other point:
        # TODO: redo for numpy arrays
        if np.all(x == np.array([0,0,0])):
            d = -np.dot(ln, oc) - np.sqrt(D)
            x = self._o + d*ln
            
        lgg.debug(x)
        
        # The vector that points from the center of the sphere to the intersection is r:
        r = x - self._c
        lgg.debug(r)
        
        # To find the angle of intersecting ray with the normal to the surface, first find the cosine of that angle (absolute value of it)
        c_angle = np.abs(np.dot(ln, r)/self._Rp)
        lgg.debug(c_angle)
        
        # And finally, return the angle (absolute value)
        if c_angle <= 1:
            pass
        elif (c_angle - 1) < 1e-5:
            c_angle = 1
        
        return np.arccos(np.abs(c_angle))
    
    # This function calculates the normalized force (i.e. actual force multiplied by c/(n_1 P)) of a single ray described by a line whose origin is o and whose direction of propagation is l. The sphere of radius R has its center in c and has refractive index nr.
    # Important note: the polarization p is a Jones' vector specified in the lab's coordinate system (e.g. before entering the lens, so that it only has XY components). This vector can be complex. For example, for circular polarization this vector would be (1,i,0), while for linear polarization it is completely real. Its normalization is not important as it is normalized in the code.
    def _ray_force(self, p):
        # Calculate the incidence angle first. If its = nan, then there is no intersection and the force is zero:
        th = self._intersection_angle()
        if np.isnan(th):
            return np.array([0,0,0])
        
        ## First we have to determine the coordinate system for the gradient and scattering forces:
        # The scattering force direction, according to Ashkin, 1992, is along the ray propagation direction:
        dir_scat = self._l/npl.norm(self._l)
        
        # The gradient force direction (Ashkin, 1992) is orthogonal to the ray propagation direction and lies in the plane formed by the ray and the center of the sphere. For that, we first make a vector that points from the center of the sphere to one of the points in the line and Gram-Schmidt orthogonalize it to make a vector perpendicular to the scattering
        a = self._o- self._c
        
        dir_grad = a - np.dot(a, dir_scat)*dir_scat
        
        # If the rays passes through the center of the sphere, then dir_grad will be = 0, but this is not a problem as this ray will not exert any gradient force. Then, we can take any direction as dir_grad without any consequence:
        if not np.all(dir_grad == np.array([0,0,0])):
            dir_grad = dir_grad / npl.norm(dir_grad)
            
        # The magnitudes of the forces are specified in Ashkin, 1992. First let's calculate some auxiliary quantities:
        
        # Refraction angle:
        r = self._snell(th)
        
        # Transmission and reflection coefficients
        # Let's calculate the projection of the polarization vector on the incidence plane and the magnitude of that projection
        Pp = (np.abs(np.dot(p, dir_grad))**2 + np.abs(np.dot(p, dir_scat))**2)/(npl.norm(p)**2)
        
        # Sometimes, the proportion will be slightly bigger than 1 because of floating-point errors. The following corrects it:
        if 1 < Pp < 1+1e-7:
            Pp = 1
        
        # Note: if dir_grad is null (when the ray is normal on the sphere), Pp will take some value between 0 and 1, but it won't matter since at normal incidence, Fresnel doesn't depend on the polarization
        
        T, R = self._fresnel(th, r, Pp)
        
        # And finally, the scattering force magnitude:
        Fs = 1 + R*np.cos(2*th) - (T**2 * (np.cos(2*th-2*r) + R*np.cos(2*th))) / (1 + R**2 + 2*R*np.cos(2*r))
        
        # And the gradient force magnitude:
        Fg = R*np.sin(2*th) - (T**2 * (np.sin(2*th-2*r) + R*np.sin(2*th))) / (1 + R**2 + 2*R*np.cos(2*r))
        
        # And calculate the total force:
        # Note that the sign of Fg is due to a sign error (or maybe misunderstanding?) in Ashkin, 1992
        F = Fs*dir_scat - Fg*dir_grad
        
        return F
    
class OpticalSystemSimpleUniform(OpticalSystem):
    def __init__(self, c, Rp, nr, Rl, f, p):
        super().__init__(c, Rp, nr)
        
        self._p = p
        
        self.set_focal_distance(f)
        self.set_lens_radius(Rl)
        
        self._c = np.array([0, 0, f]) + c
        
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
        self._c = np.array([0, 0, self._f]) + c
        
    # Returns the total force by a single ray in this case
    def _total_ray_force(self, r, th):
        self._o = np.array([r*np.cos(th), r*np.sin(th), 0])
        self._l = np.array([0, 0, self._f]) - self._o
        
        F = self._ray_force(self._p)
    
        return (r/(np.pi * self._Rl**2))*F
    
    def integrate(self):
        Ft = np.array([np.nan, np.nan, np.nan])
        
        for i in range(0,3):
            Ft[i] = si.dblquad(lambda r,th: self._total_ray_force(r, th)[i], 0, 2*np.pi, lambda x: 0, lambda x: self._Rl, epsabs=1e-4, epsrel=1e-4)[0]
        
        return Ft