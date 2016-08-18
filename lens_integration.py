import numpy as np
import scipy.integrate as si

import sphere_intersection as spi

def simple_unif(r, th, R, f):
    o = np.array([r*np.cos(th), r*np.sin(th), 0])
    l = np.array([0, 0, f]) - o
    
    return (r/(np.pi * R**2), o, l)

def simple_unif_integrate(pos, rp, n, R, f, p):
    # Force from a ray to be integrated
    def ray_force_complete(r, th):
        ray = simple_unif(r, th, R, f)
        
        P = ray[0]
        o = ray[1]
        l = ray[2]
        
        c = pos + np.array([0,0,f])
        
        F = spi.ray_force(c, rp, o, l, p, n)
        
        return F*P
    
    # Now just calculate the components of the force (warning: especially shitty code, to be refactored)
    Ft = np.array([np.nan, np.nan, np.nan])
    
    for i in range(0,3):
        Ft[i] = si.dblquad(lambda r,th: ray_force_complete(r, th)[i], 0, 2*np.pi, lambda x: 0, lambda x: R, epsabs=1e-4, epsrel=1e-4)[0]
        
    return Ft