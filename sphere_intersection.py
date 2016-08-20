import numpy as np
import numpy.linalg as npl

import logging

logging.basicConfig(level=logging.INFO)
lgg = logging.getLogger("intersection_module")

# Calculates the refraction angle given the incidence angle and the relative index of refraction
def snell(theta, nr):
    # Raise exception if the angle is out of the [0,pi/2] range
    if not 0 <= theta <= np.pi/2:
        raise ValueError("Incidence angle out of range: {0}".format(theta))
    
    return np.arcsin(1/nr*np.sin(theta))

# Calculates the transmission and reflection for a ray with a given incidence angle (th), refraction angle (r), and polarization angle (p) when the relative index of refraction is specified (nr)
# Also note that the polarization is specified as the normalized power of p-polarization (Pp). Then, the power of the s-polarization is simply (1-Pp).
def fresnel(th, r, Pp, nr):
    # Raise exception if the incidence angle is significantly greater than pi/2 or is negative
    if not 0 <= th <= np.pi/2:
        raise ValueError("Incidence angle is out of range: {0}".format(th))
    
    # The same for the transmission angle
    if not 0 <= r <= np.pi/2:
        raise ValueError("Transmission angle is out of range: {0}".format(r))
    
    # And Pp should be within the [0,1] interval
    if not 0 <= Pp <= 1:
        raise ValueError("p-polarization proportion out of range: {0}".format(Pp))
    
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
        raise ValueError("The sphere radius is not valid: {0}".format(R))
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
    
# This function calculates the normalized force (i.e. actual force multiplied by c/(n_1 P)) of a single ray described by a line whose origin is o and whose direction of propagation is l. The sphere of radius R has its center in c and has refractive index nr.
# Important note: the polarization p is a Jones' vector specified in the lab's coordinate system (e.g. before entering the lens, so that it only has XY components). This vector can be complex. For example, for circular polarization this vector would be (1,i,0), while for linear polarization it is completely real. Its normalization is not important as it is normalized in the code.
def ray_force(c, R, o, l, p, nr):
    # Calculate the incidence angle first. If its = nan, then there is no intersection and the force is zero:
    th = intersection_angle(c, R, o, l)
    if np.isnan(th):
        return np.array([0,0,0])
    
    ## First we have to determine the coordinate system for the gradient and scattering forces:
    # The scattering force direction, according to Ashkin, 1992, is along the ray propagation direction:
    dir_scat = l/npl.norm(l)
    
    # The gradient force direction (Ashkin, 1992) is orthogonal to the ray propagation direction and lies in the plane formed by the ray and the center of the sphere. For that, we first make a vector that points from the center of the sphere to one of the points in the line and Gram-Schmidt orthogonalize it to make a vector perpendicular to the scattering
    a = o-c
    
    dir_grad = a - np.dot(a, dir_scat)*dir_scat
    
    # If the rays passes through the center of the sphere, then dir_grad will be = 0, but this is not a problem as this ray will not exert any gradient force. Then, we can take any direction as dir_grad without any consequence:
    if np.all(dir_grad == np.array([0,0,0])):
       dir_grad = np.array([0,0,0])
    
    # Otherwise, we just normalize it
    else:
        dir_grad = dir_grad / npl.norm(dir_grad)
        
    # The magnitudes of the forces are specified in Ashkin, 1992. First let's calculate some auxiliary quantities:
    
    # Refraction angle:
    r = snell(th, nr)
    
    # Transmission and reflection coefficients
    # Let's calculate the projection of the polarization vector on the incidence plane and the magnitude of that projection
    Pp = (np.abs(np.dot(p, dir_grad))**2 + np.abs(np.dot(p, dir_scat))**2)/(npl.norm(p)**2)
    
    # Sometimes, the proportion will be slightly bigger than 1 because of floating-point errors. The following corrects it:
    if 1 < Pp < 1+1e-7:
        Pp = 1
    
    # Note: if dir_grad is null (when the ray is normal on the sphere), Pp will take some value between 0 and 1, but it won't matter since at normal incidence, Fresnel doesn't depend on the polarization
    
    T, R = fresnel(th, r, Pp, nr)
    
    # And finally, the scattering force magnitude:
    Fs = 1 + R*np.cos(2*th) - (T**2 * (np.cos(2*th-2*r) + R*np.cos(2*th))) / (1 + R**2 + 2*R*np.cos(2*r))
    
    # And the gradient force magnitude:
    Fg = R*np.sin(2*th) - (T**2 * (np.sin(2*th-2*r) + R*np.sin(2*th))) / (1 + R**2 + 2*R*np.cos(2*r))
    
    # And calculate the total force:
    # Note that the sign of Fg is due to a sign error (or maybe misunderstanding?) in Ashkin, 1992
    F = Fs*dir_scat - Fg*dir_grad
    
    return F