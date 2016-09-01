# This file calculates the force (adimensional factor Q) on a particle of given index, with optics of given NA in a range of x's, y's and z's.
# The calculation is done assuming that all the rays are focused in the single spot (so that there is no explicit dependence on the radius of the particle)
import optical_system as osys
import numpy as np

# Import the Python configuration file
import config

# Output file
out_file = config.out_file

# The radius of the particle will be taken as unity with no loss of generality.
Rp = config.radius
n = config.nr

# The NA will be used for calculating the lens radius, though the pair lens radius-focal length can also be specified below
# Note: if the NA is for liquid-immersion objective, then it should be divided by the index of that liquid (so that it's a number less than 1)
NA = config.NA

f = 1e5*Rp # A lot to guarantee that the particle doesn't hit the lens (nothing horrible should happen, but still)

# Finally, calculate the lens radius
Rl = f * np.tan(np.arcsin(NA))

# The first argument is the lowest x, the second is the highest x and the third is the number of steps to take.
# If not iterating over x, then set it to 0, 0, 1 (or change 0 to the desired fixed value of x)
# All the coordinates are zero when the particle is at the focus. Z decreases when the particle is closer to the lens.

xstart = config.xstart
xstop = config.xstop
xsteps = config.xsteps

if xstart == xstop:
    xsteps = 1
xs = np.linspace(xstart, xstop, xsteps)

# Same as above, for y
ystart = config.ystart
ystop = config.ystop
ysteps = config.ysteps

if ystart == ystop:
    ysteps = 1
ys = np.linspace(ystart, ystop, ysteps)

# Same as above, for z
zstart = config.zstart
zstop = config.zstop
zsteps = config.zsteps

if zstart == zstop:
    zsteps = 1
zs = np.linspace(zstart, zstop, zsteps)

# Now we generate the space of all the necessary coordinates, where every row is a position to be calculated
xx, yy, zz = np.meshgrid(xs, ys, zs)

positions = np.vstack([xx.flatten(), yy.flatten(), zz.flatten()]).transpose()

# Initialize the system (the 0,0,0 initial position is just for completeness)
opt = osys.OpticalSystemSimpleArbitrary(np.array([0,0,0]), Rp, n, Rl, f, 
                                        config.int_pol_function, **config.int_pol_arguments)

# An auxiliary function for applying on each row
def force(row):
    opt.set_particle_center(row)
    force = opt.integrate(config.rsteps, config.thsteps)
    
    return force

forces = np.apply_along_axis(force, axis=1, arr=positions)

# Make an array that includes positions and forces
out_array = np.hstack([positions, forces])
# Save it into a file
np.savetxt(out_file, out_array, delimiter="\t", fmt='%.6e')
