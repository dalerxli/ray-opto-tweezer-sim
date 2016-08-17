# Testing rig
import unittest

# Module to test
import sphere_intersection as ix

# Auxiliary
import numpy as np
import numpy.linalg as npl

# This class tests integration over all the rays coming out of a lens
class TestIntegration(unittest.TestCase):
    # Test the power coming out of lens with a uniform fill. This power has to be equal to unity (a check on the radial integration)
    def test_integrate_uniform(self):
        pass