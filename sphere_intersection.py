import numpy as np
import numpy.linalg as npl

def intersection_angle(c, R, o, l):
    # Make l (director of the line) unitary
    l = l/npl.norm(l)
    
    # Calculate the discriminant (to see whether there are any solutions)
    D = -np.dot(l, (o-c))**2 - np.dot((o-c), (o-c)) + R**2
    
    # If there are no solutions, then there is nothing else to do
    if D < 0:
        return np.nan
    
    # Otherwise, calculate the distance along the line where an intersection occurs (doesn't matter which since this is a sphere)
    d = -np.dot(l, (o-c)) + np.sqrt(D)
    
    # The point at which the intersection occurs is x:
    x = o + d*l
    
    # The vector that points from the center of the sphere to the intersection is r:
    r = x-c
    
    # To find the angle of intersecting ray with the normal to the surface, first find the cosine of that angle
    c_angle = np.dot(x, r)/(npl.norm(x)*npl.norm(r))
    
    # And finally, return the angle (absolute value)
    return np.arccos(np.abs(c_angle))
    