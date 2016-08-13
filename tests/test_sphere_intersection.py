# Testing rig
import unittest

# Module to test
import sphere_intersection as ix

# Auxiliary
import numpy as np

class SphereIntersectionTestCase(unittest.TestCase):
    def test_intersect_normal(self):
        # The function should find the intersection between a sphere (centered at c, the first argument) of radius R (second argument) with a line of origin o (third argument) and direction vector l (fourth argument)
        c = np.array([0,0,0])
        R = 1
        o = np.array([0,0,0])
        l = np.array([1,0,0])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0)
        
        # Let's now offset the sphere and the line origin
        c = np.array([1,0,0])
        o = np.array([1,0,0])
        
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
        
if __name__ == '__main__':
    unittest.main()