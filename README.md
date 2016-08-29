# ray-opto-tweezer-sim
ray-opto-tweezer-sim is a program written for simulating forces on homogeneous dielectric spheres trapped in optical tweezers where the ray optics approximation is applicable, i.e. when the particle is significantly bigger than the wavelength of the trapping light and is bigger that the waist of the light beam.

This code is written in Python3+NumPy, which should make it easy to install and use (see "Installation" below), while also being reasonably fast (courtesy of NumPy): evaluating the force on the particle at 20 particle positions per second on a single core of Intel Celeron 1007U (an entry-level laptop CPU).
