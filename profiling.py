# Modules to test
import optical_system as osys

# Auxiliary
import numpy as np
import numpy.linalg as npl

# For integration
import scipy.integrate as si 

f = 1e-3
        
# A microscope objective with NA = 1.25 (water-immersion). The half-angle of convergence is about 70 degrees
Rl = f * np.tan(np.arcsin(1.25/1.33))

# A particle of rp=5e-6. Not necessary in this case, but I'll keep for consistency.
rp = 5e-6

# And the polarization is linear
p = np.array([1,0,0])

opt = osys.OpticalSystemSimpleUniform(np.array([0,0,0]), rp, 1.5, Rl, f, p)

# Data from Ashkin, 1992
data = np.array([
    [1.2, 0.00, 0.00, 1.01*rp, -0.276, 2],
    [1.2, 0.00, 0.98*rp, 0.00, -0.313, 1],
    #[1.2, 1.05*rp, 0.00, 0.00, -0.490, 0], #This dataset is dubious: check Ashkin
    [1.4, 0.00, 0.00, 0.93*rp, -0.282, 2],
    [1.8, 0.00, 0.00, 0.88*rp, -0.171, 2]
    ])

def check(row):
    n = row[0]
    pos = row[1:4]
    targetQ = row[4]
    
    # Force index to check
    i = row[5]
    
    opt.set_particle_center(pos)
    opt.set_particle_index(n)
    
    force = opt.integrate(200, 200)
    
    return np.abs(force[int(i)] - targetQ)
    
for i in range(0,100):
    res = np.apply_along_axis(check, axis=1, arr=data)
