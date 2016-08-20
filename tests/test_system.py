# This is the transition to the OOP-based program structure (hopefully more optimized)

# Testing rig
import unittest

# Modules to test
import optical_system as osys

# Auxiliary
import numpy as np
import numpy.linalg as npl

# For integration
import scipy.integrate as si 

class OpticalSystemIntersectionTestCase(unittest.TestCase):
    def test_intersect_normal(self):
        opt = osys.OpticalSystem(np.array([0,0,0]), 1, 1.5)
        
        opt._o = np.array([0,0,0])
        opt._l = np.array([1,0,0])
        
        angle = opt._intersection_angle()
        
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_normal_offset(self):
        opt = osys.OpticalSystem(np.array([3,0,0]), 1, 1.5)
        
        opt._o = np.array([3,0,0])
        opt._l = np.array([1,0,0])
        
        angle = opt._intersection_angle()
        self.assertAlmostEqual(angle, 0)
    
    def test_intersect_normal_inv_director(self):
        opt = osys.OpticalSystem(np.array([3,0,0]), 1, 1.5)
        
        opt._o = np.array([3,0,0])
        opt._l = np.array([1,0,0])
        
        opt._l = -6 * opt._l
        
        ## Inverting the line director shouldn't have any effect:
        
        angle = opt._intersection_angle()
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_normal_angle_director(self):
        opt = osys.OpticalSystem(np.array([3,0,0]), 1, 1.5)
        
        opt._o = np.array([3,0,0])
        opt._l = np.array([1,1,1])
        
        angle = opt._intersection_angle()
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_normal_angle_director_external(self):
        # And now let the line hit the sphere normally, but at an angle to an axis
        opt = osys.OpticalSystem(np.array([5,0,0]), 1, 1.5)
        
        opt._o = np.array([0,0,5])
        opt._l = np.array([1,0,-1])
        
        angle = opt._intersection_angle()
        self.assertAlmostEqual(angle, 0, places=5)
        
    def test_intersect_normal_at_zero(self):
        # What happens if the position of the intersection is (0,0,0)?
        opt = osys.OpticalSystem(np.array([1,0,0]), 1, 1.5)
        
        opt._o = np.array([1,0,0])
        opt._l = np.array([-1,1,1])
        
        angle = opt._intersection_angle()
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_tangent(self):
        # What happens if the ray is exactly tangent to the sphere?
        opt = osys.OpticalSystem(np.array([3,0,0]), 1, 1.5)
        
        opt._o = np.array([4,0,10])
        opt._l = np.array([0,0,-1])
        
        angle = opt._intersection_angle()
        self.assertAlmostEqual(angle, np.pi/2)
        
    def test_intersect_45(self):
        # Here, the ray should intersect the sphere at 45 degrees
        opt = osys.OpticalSystem(np.array([4,0,0]), 1, 1.5)
        
        opt._o = np.array([0,0,-3])
        opt._l = np.array([1,0,1])
        
        angle = opt._intersection_angle()
        self.assertAlmostEqual(angle, np.pi/4)
        
    def test_intersect_invalid_radius(self):
        ## The system should not accept zero or negative sphere radiuses
        for R in [-1, 0]:
            with self.assertRaises(ValueError):
                opt = osys.OpticalSystem(np.array([4,0,0]), R, 1.5)
                
# Testing the helper functions to calculate the angle of refraction and the coefficients of transmission/reflection
class SystemRefractionTestCase(unittest.TestCase):
    def test_snell_normal(self):
        # At normal incidence, the angle stays zero:
        opt = osys.OpticalSystem(np.array([4,0,0]), 1, 1.5)
        r = opt._snell(0)
        
        self.assertAlmostEqual(r, 0)
        
    def test_snell_tangent(self):
        # At tangent incidence, the refr. angle should be critical:
        nr = 1.5
        opt = osys.OpticalSystem(np.array([4,0,0]), 1, nr)
        r = opt._snell(np.pi/2)
        
        self.assertAlmostEqual(np.sin(r), 1/nr)
        
    def test_snell_homogeneous(self):
        # If the relative index is 1, then there should be no change of propagation at all
        opt = osys.OpticalSystem(np.array([4,0,0]), 1, 1)
        
        for th in [0, np.pi/3, np.pi/4, np.pi/5]:
            self.assertAlmostEqual(opt._snell(th), th)
            
    def test_fresnel_normal(self):
        # Test Fresnel at normal incidence. It should not depend on the polarization, and the energy must be conserved between the transmittance and reflectance
        # Note that Fresnel only really depends on the relative index.
        # Also note that the polarization is specified as the normalized power of p-polarization (Pp). Then, the power of the s-polarization is simply (1-Pp).
        nr = 1.5
        th = 0
        r = 0
        
        opt = osys.OpticalSystem(np.array([4,0,0]), 1, nr)
        
        for Pp in [0, 0.1, 0.5, 0.8, 1]:
            T, R = opt._fresnel(th, r, Pp)
            
            self.assertAlmostEqual(T+R, 1)
            self.assertAlmostEqual(R, ((1 - nr)/(1 + nr))**2)
            
    def test_fresnel_tangent(self):
        # The reflectivity should be 1 at tangent incidence, regardless of polarization
        nr = 1.5
        th = np.pi/2
        
        opt = osys.OpticalSystem(np.array([4,0,0]), 1, nr)
        
        r = opt._snell(th)
        
        for Pp in [0, 0.1, 0.5, 0.8, 1]:
            T, R = opt._fresnel(th, r, Pp)
            
            self.assertAlmostEqual(T+R, 1)
            self.assertAlmostEqual(R, 1)
            
    def test_fresnel_brewster(self):
        # The reflectivity of a p-polarized ray should be nearly 0 at Brewster's angle
        nr = 1.5
        th = np.arctan(nr)
        
        opt = osys.OpticalSystem(np.array([4,0,0]), 1, nr)
        
        r = opt._snell(th)
        
        T, R = opt._fresnel(th, r, 1)
        self.assertAlmostEqual(R, 0)
        
# Test of a higher-level function that returns the Q-force (explained below) of a single ray incident on a sphere
# The Q force is the force exerted by a single ray, but divided by (P*n_1/c) to not depend on the power of this ray and other factors that will have to be included later (see Ashkin, 1992)
class SphereIntersectionForceTestCase(unittest.TestCase):
    def test_force_ashkin(self):
        # Test the maximum gradient forces published in Ashkin, 1992
        
        # The origin of the ray will be on the sphere for simpler calculations:
        c = np.array([1,0,0])
        R = 1
        o = np.array([0,0,0])
        
        opt = osys.OpticalSystem(np.array([4,0,0]), 1, 1.5)
        opt._c = c
        
        # And this is the polarization of the ray in Jones notation (circular)
        p = np.array([1,1j,0])
        
        data = np.array([
            [1.1, np.sqrt(0.429**2 + 0.262**2), 79*np.pi/180],
            [1.2, np.sqrt(0.506**2 + 0.341**2), 72*np.pi/180],
            [1.4, np.sqrt(0.566**2 + 0.448**2), 64*np.pi/180],
            [1.6, np.sqrt(0.570**2 + 0.535**2), 60*np.pi/180],
            [1.8, np.sqrt(0.547**2 + 0.625**2), 59*np.pi/180],
            [2.0, np.sqrt(0.510**2 + 0.698**2), 59*np.pi/180],
            [2.5, np.sqrt(0.405**2 + 0.837**2), 64*np.pi/180]
            ])
        
        def check(row):
            th = row[2]
            Q = row[1]
            nr = row[0]
            
            opt.set_particle_index(nr)
            
            opt._l = np.array([np.cos(th), 0, np.sin(th)])
            force = opt._ray_force(p)
            
            return np.abs(npl.norm(force) - Q)
            
        res = np.apply_along_axis(check, axis=1, arr=data)
        self.assertLess(np.max(res), 0.021)
        
# This class tests integration over all the rays coming out of a lens
class TestIntegration(unittest.TestCase):
    # Finally, compare the total forces on the sphere to some of Ashkin's results
    def test_force_uniform(self):
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
            #[1.2, 1.05*rp, 0.00, 0.00, -0.490, 0],
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
            
            force = opt.integrate()
            
            return np.abs(force[int(i)] - targetQ)
            
        res = np.apply_along_axis(check, axis=1, arr=data)
        self.assertLess(np.max(res), 0.01)