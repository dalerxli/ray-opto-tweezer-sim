# Testing rig
import unittest

# Module to test
import sphere_intersection as ix

# Auxiliary
import numpy as np

# The function should find the intersection between a sphere (centered at c, the first argument) of radius R (second argument) with a line of origin o (third argument) and direction vector l (fourth argument)
class SphereIntersectionAngleTestCase(unittest.TestCase):
    def test_intersect_normal(self):
        c = np.array([0,0,0])
        R = 1
        o = np.array([0,0,0])
        l = np.array([1,0,0])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_normal_offset(self):
        # Let's now offset the sphere and the line origin
        c = np.array([3,0,0])
        o = np.array([3,0,0])
        R = 1
        l = np.array([1,0,0])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0)
    
    def test_intersect_normal_inv_director(self):
        c = np.array([3,0,0])
        o = np.array([3,0,0])
        R = 1
        l = np.array([1,0,0])
        
        # Inverting the line director shouldn't have any effect:
        l = -6*l
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_normal_angle_director(self):
        # A director can go in any direction from the center of the sphere without changing the incidence angle
        c = np.array([3,0,0])
        o = np.array([3,0,0])
        R = 1
        l = np.array([1,1,1])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_normal_angle_director_external(self):
        # And now let the line hit the sphere normally, but at an angle to an axis
        c = np.array([5,0,0])
        o = np.array([0,0,5])
        l = np.array([1,0,-1])
        R = 1
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0, places=5)
        
    def test_intersect_normal_at_zero(self):
        # What happens if the position of the intersection is (0,0,0)?
        c = np.array([1,0,0])
        R = 1
        o = np.array([1,0,0])
        l = np.array([-1,0,0])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_tangent(self):
        # What happens if the ray is exactly tangent to the sphere?
        c = np.array([3,0,0])
        R = 1
        o = np.array([4,0,10])
        l = np.array([0,0,-1])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, np.pi/2)
        
    def test_intersect_45(self):
        # Here, the ray should intersect the sphere at 45 degrees
        c = np.array([4,0,0])
        R = 1
        o = np.array([0,0,-3])
        l = np.array([1,0,1])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, np.pi/4)
        
    def test_intersect_invalid_radius(self):
        # The function should not accept zero or negative sphere radiuses
        c = np.array([4,0,0])
        o = np.array([0,0,-3])
        l = np.array([1,0,1])
        
        for R in [-1, 0]:
            with self.assertRaises(ValueError):
                ix.intersection_angle(c, R, o, l)

# Testing the helper functions to calculate the angle of refraction and the coefficients of transmission/reflection
class RefractionTestCase(unittest.TestCase):
    def test_snell_normal(self):
        # At normal incidence, the angle stays zero:
        nr = 1.5
        r = ix.snell(0, nr)
        
        self.assertAlmostEqual(r, 0)
        
    def test_snell_tangent(self):
        # At tangent incidence, the refr. angle should be critical:
        nr = 1.5
        r = ix.snell(np.pi/2, nr)
        
        self.assertAlmostEqual(np.sin(r), 1/nr)
        
    def test_snell_homogeneous(self):
        # If the relative index is 1, then there should be no change of propagation at all
        nr = 1
        
        for th in [0, np.pi/3, np.pi/4, np.pi/5]:
            self.assertAlmostEqual(ix.snell(th, nr), th)
            
    def test_snell_invalid_angles(self):
        # Snell shouldn't accept any negative angles or angles greater than pi/2
        nr = 1.5
        
        for th in [np.pi/1.5, -np.pi/3]:
            with self.assertRaises(ValueError):
                ix.snell(th, nr)
            
    def test_fresnel_normal(self):
        # Test Fresnel at normal incidence. It should not depend on the polarization, and the energy must be conserved between the transmittance and reflectance
        # Note that Fresnel only really depends on the relative index.
        # Also note that the polarization is specified as the normalized power of p-polarization (Pp). Then, the power of the s-polarization is simply (1-Pp).
        nr = 1.5
        th = 0
        r = 0
        
        for Pp in [0, 0.1, 0.5, 0.8, 1]:
            T, R = ix.fresnel(th, r, Pp, nr)
            
            self.assertAlmostEqual(T+R, 1)
            self.assertAlmostEqual(R, ((1 - nr)/(1 + nr))**2)
            
    def test_fresnel_tangent(self):
        # The reflectivity should be 1 at tangent incidence, regardless of polarization
        nr = 1.5
        th = np.pi/2
        r = ix.snell(th, nr)
        
        for Pp in [0, 0.1, 0.5, 0.8, 1]:
            T, R = ix.fresnel(th, r, Pp, nr)
            
            self.assertAlmostEqual(T+R, 1)
            self.assertAlmostEqual(R, 1)
            
    def test_fresnel_brewster(self):
        # The reflectivity of a p-polarized ray should be nearly 0 at Brewster's angle
        nr = 1.5
        th = np.arctan(nr)
        
        r = ix.snell(th, nr)
        
        T, R = ix.fresnel(th, r, 1, nr)
        self.assertAlmostEqual(R, 0)
        
    def test_fresnel_invalid_arguments(self):
        nr = 1.5
        th = np.pi/4
        r = ix.snell(th, nr)
        
        # Fresnel shouldn't accept Pp out of [0,1] interval
        for Pp in [-0.1, 1.2]:
            with self.assertRaises(ValueError):
                ix.fresnel(th, r, Pp, nr)
            
        # Or angles that are too big or negative
        for th in [-0.1, np.pi/1.5]:
            for r in [-0.1, np.pi/1.5]:
                with self.assertRaises(ValueError):
                    ix.fresnel(th, r, 0.5, nr)
        
# Test of a higher-level function that returns the Q-force (explained below) of a single ray incident on a sphere
# The Q force is the force exerted by a single ray, but divided by (P*n_1/c) to not depend on the power of this ray and other factors that will have to be included later (see Ashkin, 1992)
class SphereIntersectionForceTestCase(unittest.TestCase):
    def test_force_normal(self):
        # The ray incides normally on the sphere
        c = np.array([5,0,0])
        R = 1
        o = np.array([0,0,0])
        l = np.array([1,0,0])
        
        # This is the relative refractive index of the sphere
        nr = 1.5
        
        force = ix.ray_force(c, R, o, l, nr)
        
        # The force should be in the direction of the ray in this case
        self.assertGreater(force[0], 0)
        
        # And it should not have any orthogonal components
        self.assertGreater(force[1], 0)
        self.assertGreater(force[2], 0)
        
    def test_force_no_particle(self):
        # If the particle has relative index = 1, then there is no particle and no force
        c = np.array([5,0,0])
        R = 1
        o = np.array([0,0,0])
        l = np.array([1,0,0])
        
        # This is the relative refractive index of the sphere
        nr = 1
        
        force = ix.ray_force(c, R, o, l, nr)
        
        # There should be no force
        self.assertAlmostEqual(np.norm(force), 0)
        
if __name__ == '__main__':
    unittest.main()