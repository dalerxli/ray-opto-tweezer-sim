# Testing rig
import unittest

# Module to test
import sphere_intersection as ix

# Auxiliary
import numpy as np

# The function should find the intersection between a sphere (centered at c, the first argument) of radius R (second argument) with a line of origin o (third argument) and direction vector l (fourth argument)
class SphereIntersectionTestCase(unittest.TestCase):
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
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_normal_at_zero(self):
        # What happens if the position of the intersection is (0,0,0)?
        c = np.array([1,0,0])
        R = 1
        o = np.array([1,0,0])
        l = np.array([-1,0,0])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0)
        
    def test_intersect_tangent(self):
        c = np.array([3,0,0])
        R = 1
        o = np.array([4,0,10])
        l = np.array([0,0,-1])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, np.pi/2)
        
if __name__ == '__main__':
    unittest.main()