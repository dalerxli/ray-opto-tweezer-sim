import numpy as np
import numpy.linalg as npl

import logging

logging.basicConfig(level=logging.INFO)
lgg = logging.getLogger("intersection_module")

# Calculates the refraction angle given the incidence angle and the relative index of refraction
def snell(theta, nr):
    # Raise exception if the angle is out of the [0,pi/2] range
    if not 0 <= theta <= np.pi/2:
        raise ValueError("Incidence angle out of range")
    
    return np.arcsin(1/nr*np.sin(theta))

# Calculates the transmission and reflection for a ray with a given incidence angle (th), refraction angle (r), and polarization angle (p) when the relative index of refraction is specified (nr)
# Also note that the polarization is specified as the normalized power of p-polarization (Pp). Then, the power of the s-polarization is simply (1-Pp).
def fresnel(th, r, Pp, nr):
    # Raise exception if the incidence angle is significantly greater than pi/2 or is negative
    if not 0 <= th <= np.pi/2:
        raise ValueError("Incidence angle is out of range")
    
    # The same for the transmission angle
    if not 0 <= r <= np.pi/2:
        raise ValueError("Transmission angle is out of range")
    
    # And Pp should be within the [0,1] interval
    if not 0 <= Pp <= 1:
        raise ValueError("p-polarization proportion out of range")
    
    # Calculate the reflectivities:
    Rs = ((np.cos(th) - nr*np.cos(r))/(np.cos(th) + nr*np.cos(r)))**2
    Rp = ((np.cos(r) - nr*np.cos(th))/(np.cos(r) + nr*np.cos(th)))**2
    
    # Calculate the final reflectivity and transmittivity
    R = Rs*(1-Pp) + Rp*Pp
    T = 1 - R
    
    # And return the value
    return (T, R)

def intersection_angle(c, R, o, l):
    # Make sure that the sphere radius is not zero or negative
    if R <= 0:
        raise ValueError("The sphere radius is negative or zero")
    # Make l (director of the line) unitary
    ln = l/npl.norm(l)
    
    oc = (o-c)
    
    # Calculate the discriminant (to see whether there are any solutions)
    D = np.dot(ln, oc)**2 - np.dot(oc, oc) + R**2
    lgg.debug(D)
    
    # If there are no solutions, then there is nothing else to do
    if D < 0:
        return np.nan
    
    # Otherwise, calculate the distance along the line where an intersection occurs (doesn't matter which since this is a sphere)
    d = -np.dot(ln, oc) + np.sqrt(D)
    lgg.debug(d)
    
    # The point at which the intersection occurs is x:
    x = o + d*ln
    
    # If x is zero (which would be a problem when calculating its inverse norm), switch to the other point:
    if np.all(x == np.array([0,0,0])):
        d = -np.dot(ln, oc) - np.sqrt(D)
        x = o + d*ln
        
    lgg.debug(x)
    
    # The vector that points from the center of the sphere to the intersection is r:
    r = x-c
    lgg.debug(r)
    
    # To find the angle of intersecting ray with the normal to the surface, first find the cosine of that angle (absolute value of it)
    c_angle = np.abs(np.dot(ln, r)/R)
    lgg.debug(c_angle)
    
    # And finally, return the angle (absolute value)
    if c_angle <= 1:
        pass
    elif (c_angle - 1) < 1e-5:
        c_angle = 1
    return np.arccos(np.abs(c_angle))
    
    