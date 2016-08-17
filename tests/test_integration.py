# Testing rig
import unittest

# Modules to test
import sphere_intersection as ix
import lens_integration as lint

# Auxiliary
import numpy as np
import numpy.linalg as npl

# For integration
import scipy.integrate as si

# This class tests integration over all the rays coming out of a lens
class TestIntegration(unittest.TestCase):
    # Test the power coming out of lens with a uniform fill. This power has to be equal to unity (a check on the polar integration)
    def test_integrate_simple_uniform(self):
        # Radius of the lens
        R = 1e-3
        
        # Its focal distance (just for completeness)
        f = 1e-3
        
        # A helper function to call the function from lens_integration since that function requires more arguments
        def unif_inten(r, th):
            # The target function should return the intensity multiplied by r (for polar integration) as the first element
            return lint.simple_unif(r, th, R, f)[0]
        
        int_result = si.dblquad(unif_inten, 0, 2*np.pi, lambda x: 0, lambda x: R)
        
        self.assertLess(int_result[1], 1e-4)
        self.assertAlmostEqual(int_result[0], 1)
        
    # Test the direction of rays that come out of the lens with the simple uniform fill. They should intersect at a single point
    def test_direction_simple_uniform(self):
        # Radius of the lens
        R = 1e-3
        
        # Its focal distance
        f = 1e-3
        
        # Now check some rays coming out of the lens
        r_space = np.linspace(0, R, 20)
        th_space = np.linspace(0, 2*np.pi, 20)
        
        # Note: speed is not too important here, so the following is fine)
        for r in r_space:
            for th in th_space:
                ray = lint.simple_unif(r, th, R, f)
                # The second returned value should be the origin of the ray
                o = ray[1]
                
                self.assertFalse(np.any(np.isnan(o)))
                
                # And the third is the direction
                l = ray[2]
                self.assertFalse(np.any(np.isnan(l)))
                
                # Normalize it just in case
                l = l/npl.norm(l)
                
                # At a distance of sqrt(r**2 + f**2) there should be an intersection at (0,0,f), assuming that the center of the lens is at (0,0,0) and the plane of the lens is perpendicular to z
                d = np.sqrt(r**2 + f**2)
                
                self.assertLess(npl.norm(o + d*l - np.array([0,0,f])), 1e-8)
    
if __name__ == '__main__':
    unittest.main()