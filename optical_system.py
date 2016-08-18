import numpy as np
import numpy.linalg as npl

import logging

logging.basicConfig(level=logging.INFO)
lgg = logging.getLogger("intersection_module")

class OpticalSystem(object):
    def __init__(self, c, Rp, nr):
        # Particle properties
        self._nr = nr
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
    
    # Calculates the refraction angle given the incidence angle and the relative index of refraction
    def _snell(self, theta):
        # Raise exception if the angle is out of the [0,pi/2] range
        if not 0 <= theta <= np.pi/2:
            raise ValueError("Incidence angle out of range: {0}".format(theta))
        
        return np.arcsin(1/self._nr*np.sin(theta))
    
    def _intersection_angle(self):
        # Make l (director of the line) unitary
        ln = self._l/npl.norm(self._l)
        
        oc = self._o - self._c
        
        # Calculate the discriminant (to see whether there are any solutions)
        D = np.dot(ln, oc)**2 - np.dot(oc, oc) + self._Rp**2
        lgg.debug(D)
        
        # If there are no solutions, then there is nothing else to do
        if D < 0:
            return np.nan
        
        # Otherwise, calculate the distance along the line where an intersection occurs (doesn't matter which since this is a sphere)
        d = -np.dot(ln, oc) + np.sqrt(D)
        lgg.debug(d)
        
        # The point at which the intersection occurs is x:
        x = self._o + d*ln
        
        # If x is zero (which would be a problem when calculating its inverse norm), switch to the other point:
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