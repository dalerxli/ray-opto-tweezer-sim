# Testing rig
import unittest

# Module to test
import sphere_intersection as ix

# Auxiliary
import numpy as np

class SphereIntersectionTestCase(unittest.TestCase):
    def test_intersect_zero(self):
        # The function should find the intersection between a sphere (centered at c, the first argument) of radius R (second argument) with a line of origin o (third argument) and direction vector l (fourth argument)
        c = np.array([0,0,0])
        R = 1
        o = np.array([0,0,0])
        l = np.array([1,0,0])
        
        angle = ix.intersection_angle(c, R, o, l)
        self.assertAlmostEqual(angle, 0)
        
if __name__ == '__main__':
    unittest.main()