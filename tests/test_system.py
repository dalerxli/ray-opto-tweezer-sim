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