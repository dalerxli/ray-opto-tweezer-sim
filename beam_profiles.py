# Import the necessary packages for the correct calculations
from numpy import pi
import numpy as np

### IMPORTANT NOTE: in your functions, always use the functions provided by numpy or compatible packages. For that, use e.g. np.sin(x), np.exp(x) etc.

# This one is Gaussian with fixed polarization (spatially uniform polarization). Note how we access the optional arguments for the function
# This function requires a 'p' argument that is a 2D numpy array that specifies the polarization in Jones' notation
def gaussian_fixed(r, th, Rl, **kwargs):
    # The total number of rays in the simulation is necessary for constructing the right arrays later. Don't worry about it
    n_rays = len(r)
    
    # The polarization is an external argument
    p = np.hstack([kwargs['p'], 0])
    
    # So is the beam width at the lens, though we multiply it by the radius of the lens (easier to specify)
    a = kwargs['a']*Rl
    
    # First, we calculate the normalization:
    # The total power passing through the aperture of the lens must be 1. Then, the total power of the full beam is
    P_0 = 1/ (1 - np.exp(-2 * (Rl/a)**2) )
    
    # The peak intensity is then:
    I_0 = 2*P_0/(pi * a**2)
    
    # And the normalized intensity function is
    I = I_0 * np.exp(-2 * (r/a)**2)
    
    # The polarization is fixed. We just make an array of appropriate dimensions here, no need to modify it
    pol = np.tile(p, (n_rays, 1))
    
    return np.hstack([I.reshape(-1,1), pol])

# This one is Gaussian with radial polarization
def gaussian_radial(r, th, Rl, **kwargs):
    # The total number of rays in the simulation is necessary for constructing the right arrays later. Don't worry about it
    n_rays = len(r)
    
    # So is the beam width at the lens
    a = kwargs['a']
    
    # First, we calculate the normalization:
    # The total power passing through the aperture of the lens must be 1. Then, the total power of the full beam is
    P_0 = 1/ (1 - np.exp(-2 * (Rl/a)**2) )
    
    # The peak intensity is then:
    I_0 = 2*P_0/(pi * a**2)
    
    # And the normalized intensity function is
    I = I_0 * np.exp(-2 * (r/a)**2)
    
    # The polarization is radial. We just make an array of appropriate dimensions here, no need to modify it
    pol = np.array([np.cos(th), np.sin(th), np.zeros(r.shape)]).transpose()
    
    return np.hstack([I.reshape(-1,1), pol])