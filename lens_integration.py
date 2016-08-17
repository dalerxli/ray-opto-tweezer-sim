import numpy as np

def simple_unif(r, th, R, f):
    o = np.array([r*np.cos(th), r*np.sin(th), 0])
    l = np.array([0, 0, f]) - o
    
    return (r/(np.pi * R**2), o, l)
