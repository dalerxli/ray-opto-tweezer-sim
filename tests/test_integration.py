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
    def test_integrate_uniform(self):
        # Radius of the lens
        R = 1e-3
        
        # A helper function to call the function from lens_integration since that function requires more arguments
        def unif_inten(r, th):
            # The target function should return the intensity multiplied by r (for polar integration) as the first element
            return lint.unif(r, th, R)[0]
        
        int_result = si.dblquad(unif_inten, 0, 2*np.pi, lambda x: 0, lambda x: R)
        
        self.assertLess(int_result[1], 1e-4)
        self.assertAlmostEqual(int_result[0], 1)
        
if __name__ == '__main__':
    unittest.main()