import beam_profiles as bp
# This file serves to specify the parameters of the simulation

# Miscellaneous options
# File into which the computed data will be written
out_file = "results.tsv"

### Particle settings
# The radius of the particle. When the rays are focused in a single spot, it's not important (and can be set to unity for easier data interpretation), but for aberrated beams this does have an effect.
radius = 1

# The relative refractive index of the particle (index of the particle divided by the index of the medium)
nr = 1.5

### Beam settings

# First, we select an intensity/polarization function
int_pol_function = bp.gaussian_fixed

# And then specify the optional arguments for that function (like beam waist and polarization)
int_pol_arguments = {'a': 1.0, 'p': np.array([1,0])}

# Numerical aperture of the lens that is used to focus light onto the particle
NA = 0.85

### Integration settings
# Number of steps into which the radial coordinate of the lens will be subdivided for integration
rsteps = 200

# Number of steps into which the azimuthal coordinate of the lens will be subdivided for integration
thsteps = 200

### Position settings
# The range of positions (for each coordinate) on which the force will be calculated. The positions are relative to the focal point, and negative Z is closer to the lens. The positions are dimensional (i.e. measured in meters or whichever units you are using). It can be handy to set the particle radius to unity in order to have the positions in terms of it (which can be done without losing generality when all the rays are focused into a single spot).

# If you wouldn't like to vary the x coordinate, just set xstart and xstop to the same value. xsteps will then be ignored
xstart = 0
xstop = 0
xsteps = 100

# Absolutely analogous to the above
ystart = 0
ystop = 0
ysteps = 100

# Absolutely analogous to the above
zstart = -2
zstop = 2
zsteps = 100
