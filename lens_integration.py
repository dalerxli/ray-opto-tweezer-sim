import numpy as np
import sphere_intersection as spi

def simple_unif(r, th, R, f):
    o = np.array([r*np.cos(th), r*np.sin(th), 0])
    l = np.array([0, 0, f]) - o
    
    return (r/(np.pi * R**2), o, l)

def simple_unif_integrate(x, y, z, rp, n, R, f, p):
    pass