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
        
# Test of a higher-level function that returns the Q-force (explained below) of a single ray incident on a sphere
# The Q force is the force exerted by a single ray, but divided by (P*n_1/c) to not depend on the power of this ray and other factors that will have to be included later (see Ashkin, 1992)
#class SphereIntersectionForceTestCase(unittest.TestCase):
    #def test_force_normal(self):
        ## The ray incides normally on the sphere
        #c = np.array([5,0,0])
        #R = 1
        #o = np.array([0,0,0])
        #l = np.array([1,0,0])
        
        ## This is the absolute refractive index of the sphere
        #n = 1.5
        
        #force = ix.ray_force(c, R, o, l, n)
        
        ## The force should be in the direction of the ray in this case
        #self.assertGreater(force[0], 0)
        
        ## And it should not have any orthogonal components
        #self.assertGreater(force[1], 0)
        #self.assertGreater(force[2], 0)
if __name__ == '__main__':
    unittest.main()