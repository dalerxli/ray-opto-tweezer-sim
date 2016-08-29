# ray-opto-tweezer-sim
##Short description
`ray-opto-tweezer-sim` is a program written for simulating forces on homogeneous dielectric spheres trapped in optical tweezers where the ray optics approximation is applicable, i.e. when the particle is significantly bigger than the wavelength of the trapping light and is bigger that the waist of the light beam.

This code is written in Python3+NumPy, which should make it easy to install and use (see "Installation" below), while also being reasonably fast (courtesy of NumPy): evaluating the force on a particle is done at 20 particle positions per second on a single core of Intel Celeron 1007U (an entry-level laptop CPU).

The theoretical background for this code is described in Ashkin, 1992. Since I couldn't find the code used by Ashkin in his paper, I have decided to write this one.

### Current features
- Calculation of the Q-factor (see Ashkin, 1992 for the definition) is done at a relatively high speed and with controllable precision. The Q factor is computed relatively to the power that passes through the lens (not the total incident power on the lens), as this approach seems to be more convenient.

- The program (now) accepts a TEM00 Gaussian beam with an arbitrary beam waist specified by the ratio between the beam waist and the lens radius. That way, by specifying the NA of the focusing lens, the cone of light is completely determined. In case if the Gaussian beam is significantly larger than the lens, it is clipped. That way, a uniform lens illumination can be considered by setting a to a large value (say, 1000), so that uniform illumination cases are also considered.

- Any spatially uniform polarization state of the incoming light is possible, such as linear polarization (e.g. *(1,0)* in Jones notation), circular (equivalent here to non-polarized), specified by *(1,1j)* where *j* is the imaginary unit. More general elliptic polarization states are also possible (equivalent here to arbitrary mixtures of p and s polarizations).

- Easy setting-up: all the simulation values and the output file are set in the amply-commented configuration file "config.ini".

- Particle positions to be evaluated are specified as ranges in each of the X,Y,Z coordinates, which allows to construct 1D, 2D and 3D force profiles with an arbitrary number of steps in each of the coordinates, while increasing the performance (some calculations can be recycled by the program).

- Trustworthy calculations: an automated test suite (in `tests` subdirectory) verifies the operation of the code, beginning with the basics (Snell law and Fresnel coefficients implementation) and ending by the total force exerted on a particle. Where possible, the test suite ensures that the generated results are consistent with the published ones (again, with Ashkin, 1992). Special attention is paid to corner cases for each of the simulation functions in order to guarantee their correctness.

### Upcoming features
- Arbitrary beam intensity profile that can be easily specified by the user.

- Arbitrary spatial polarization profile to allow simulating e.g. radial polarization.

- Multithreading support to speed things up on multi-core or multi-CPU machines.

- Arbitrary ray-transfer function for the lens to allow simulating e.g. lenses with aberrations or lenses more interesting than the occasional Abbe-compliant objectives.